#include "unit_tester.hpp"

#include "graphics/new/vao.h"
#include "graphics/new/window.h"

#include "graphics/colormap.h"
#include "graphics/shader.h"

#include "library/ckf.h"

#include "mesh/field.hpp"
#include "mesh/points.h"
#include "mesh/topology.h"

using namespace avro;
using namespace avro::graphics;

UT_TEST_SUITE( graphics_decomposition_suite )

class TestField : public Field<Simplex,std::vector<real_t>> {

public:
  TestField( const Topology<Simplex>& topology , coord_t order ) :
    Field( topology , order, CONTINUOUS )
  {
    build();

    const Simplex& element = this->element();

    index_t n = 0;
    for (index_t k = 0; k < topology.nb(); k++) {
      for (index_t j = 0; j < element.nb_basis(); j++) {

        // retrieve the reference coordinate
        const real_t* xref = element.reference().get_reference_coordinate(j);

        // evaluate the basis functions of the mesh element
        std::vector<real_t> phi( topology.element().nb_basis() );
        topology.element().reference().basis().evaluate( xref , phi.data() );

        std::vector<real_t> x(4,0);
        for (index_t i = 0; i < phi.size(); i++) {
          x[0] += phi[i]*topology.points()[ topology(k,i) ][0];
          x[1] += phi[i]*topology.points()[ topology(k,i) ][1];

          if (topology.points().dim() > 2)
            x[2] += phi[i]*topology.points()[ topology(k,i) ][2];
        }

        index_t idx = this->index(k,j);
        if (this->type() == DISCONTINUOUS) {
          avro_assert( n++ == idx );
        }

        x[3] = 0.95*sin( x[0]*M_PI )*sin( x[1]*M_PI )* cos( x[2]*M_PI );
        this->value(idx) = x;
      }
    }
  }

  index_t nb_rank() const { return 4; }
};

graphics::mat4 model_matrix = glm::identity();

index_t width = 400;
index_t height = 400;
bool dragging = false;
int xm = -1;
int ym = -1;
int XM,YM;
VertexAttributeObject* vao_ptr = nullptr;
ShaderProgram* triangle_shader = nullptr;
ShaderProgram* edge_shader = nullptr;
ShaderProgram* point_shader = nullptr;
GLFWwindow* win = nullptr;

// set up the view
graphics::vec3 eye = {0,0,5};
graphics::vec3 up = {0,1,0};
graphics::vec3 center = {0.5,0.5,0.5};

float fov    = M_PI/4.0;
float aspect = width/height;
float znear  = 0.001;
float zfar   = 100;

// define the transformation matrices, and set the MVP into the shader program
graphics::mat4 perspective_matrix = glm::perspective( fov , aspect , znear , zfar );
graphics::mat4 view_matrix = glm::lookAt( eye , center , up );
graphics::mat4 mv = view_matrix * model_matrix;
graphics::mat4 mvp = perspective_matrix * mv;

mat4
rotation( real_t X , real_t Y ) {
  real_t X2 = X*X, Y2 = Y*Y;
  real_t q = 1 + X2 + Y2;
  real_t s = 1 - X2 - Y2;
  real_t r2 = 1/(q*q), s2 = s*s;
  real_t A = (s2 + 4*(Y2 - X2))*r2;
  real_t B = -8*X*Y*r2;
  real_t C = 4*s*X*r2;
  real_t D = (s2 + 4*(X2 - Y2))*r2;
  real_t E = 4*s*Y*r2;
  real_t F = (s2 - 4*(X2 + Y2))*r2;

  mat4 R; // initializes to zero
  R(0,0) =  A; R(1,0) =  B; R(2,0) = C;
  R(0,1) =  B; R(1,1) =  D; R(2,1) = E;
  R(0,2) = -C; R(1,2) = -E; R(2,2) = F;
  R(3,3) =  1;
  return R;
}

// for now, hack the center of rotation to be the center of the cube domain
real_t xc = 0.5;
real_t yc = 0.5;
real_t zc = 0.5;
mat4 T = glm::translate( glm::identity() , {xc,yc,zc} );
mat4 Tinv = glm::translate( glm::identity() , {-xc,-yc,-zc} );

void
draw() {
  float col = 1.0;
  glClearColor(col,col,col, 0.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  mv  = view_matrix * model_matrix;
  mvp = perspective_matrix * mv;

  triangle_shader->use();
  triangle_shader->setUniform("u_ModelViewProjectionMatrix",mvp);

  edge_shader->use();
  edge_shader->setUniform("u_ModelViewProjectionMatrix",mvp);

  point_shader->use();
  point_shader->setUniform("u_ModelViewProjectionMatrix",mvp);

  vao_ptr->draw_edges(*edge_shader);
  vao_ptr->draw_triangles(*triangle_shader);
  //vao_ptr->draw_points(*point_shader);

  glfwSwapBuffers(win);

}

static void
mouse_button_callback(GLFWwindow* window , int button, int action, int mods)
{
  if (action == GLFW_PRESS) dragging = true;
  else {
    xm = ym = -1;
    dragging = false;
  }
}

static void
mouse_move_callback(GLFWwindow* window, double x, double y)
{
  if (!dragging) { return; }
  if (xm < 0 || ym < 0) {
    xm = x;
    ym = y;
    return;
  }

  // translate to the center, rotate, then translate back
  mat4 R = rotation( -(x - xm)/width , (y - ym)/height );
  model_matrix = T * R * Tinv * model_matrix;

  draw();

  xm = x;
  ym = y;
}

void
change_rank(index_t rank) {
  vao_ptr->set_rank(rank);
  draw();
}

UT_TEST_CASE( simplices_2d_test )
{
  coord_t number = 3;
  coord_t dim = number;
  std::vector<index_t> dims(number,10);
  CKF_Triangulation topology( dims );

  coord_t geometry_order = 2;
  Points nodes(dim);
  topology.element().set_basis( BasisFunctionCategory_Lagrange );
  Topology<Simplex> curvilinear(nodes,topology,geometry_order);
  curvilinear.element().set_basis( BasisFunctionCategory_Lagrange );

  #if 1
  for (index_t k = 0; k < curvilinear.points().nb(); k++) {

    real_t s = curvilinear.points()[k][0];
    real_t t = curvilinear.points()[k][1];

    real_t theta = s*M_PI;
    real_t R = 0.75 + t*(1.0 - 0.75);

    curvilinear.points()[k][0] = R*cos(theta);
    curvilinear.points()[k][1] = R*sin(theta);
  }
  #endif

  coord_t solution_order = 3;
  std::shared_ptr<TestField> field = std::make_shared<TestField>(curvilinear,solution_order);
  field->element().set_basis( BasisFunctionCategory_Lagrange );
  curvilinear.fields().make( "test" , field );

  Window window(width,height);
  window.init();
  win = window.window();

  VertexAttributeObject vao(number,1);
  vao.build(curvilinear);

  std::vector<std::string> macros = {"#define SOLUTION_ORDER " + std::to_string(field->element().order()),
                                     "#define GEOMETRY_ORDER " + std::to_string(curvilinear.element().order()) };

  bool with_tess = true;
  ShaderProgram tshader("triangles",with_tess,macros);
  tshader.use();
  tshader.setUniform("u_ModelViewProjectionMatrix",mvp);

  ShaderProgram eshader("edges",with_tess,macros);
  eshader.use();
  eshader.setUniform("u_ModelViewProjectionMatrix",mvp);

  ShaderProgram pshader("points");
  pshader.use();
  pshader.setUniform("u_ModelViewProjectionMatrix",mvp);

  glfwSetCursorPosCallback(window.window(),&mouse_move_callback);
  glfwSetMouseButtonCallback(window.window(),&mouse_button_callback);

  // bind the colormap values to a buffer
  gl_index colormap_buffer;
  GL_CALL( glGenBuffers( 1 , &colormap_buffer ) );
  Colormap colormap;
  colormap.change_style("viridis");
  index_t ncolor = 256*3;
  GL_CALL( glBindBuffer( GL_TEXTURE_BUFFER , colormap_buffer) );
  GL_CALL( glBufferData( GL_TEXTURE_BUFFER , sizeof(float) * ncolor , colormap.data() , GL_STATIC_DRAW) );

  // generate a texture to hold the colormap buffer
  GLuint colormap_texture;
  GL_CALL( glGenTextures( 1 , &colormap_texture) );
  GL_CALL( glActiveTexture( GL_TEXTURE0 + 1 ) );
  GL_CALL( glBindTexture( GL_TEXTURE_BUFFER , colormap_texture) );
  GL_CALL( glTexBuffer( GL_TEXTURE_BUFFER , GL_R32F , colormap_buffer ) );

  // set the tessellation level for the TCS
  int level = 4;

  tshader.use();
  tshader.setUniform( "u_level" , level );

  eshader.use();
  eshader.setUniform("u_level" , level );

  vao_ptr = &vao;
  triangle_shader = &tshader;
  edge_shader = &eshader;
  point_shader = &pshader;
  draw();

  while (true) {

    // wait for user input
    glfwWaitEvents();

    // determine if we should exit the render loop
    if (glfwWindowShouldClose(window.window())) break;
    if (glfwGetKey(window.window(), GLFW_KEY_ESCAPE ) == GLFW_PRESS) break;

    if (glfwGetKey(window.window() , GLFW_KEY_0) == GLFW_PRESS) change_rank(0);
    if (glfwGetKey(window.window() , GLFW_KEY_1) == GLFW_PRESS) change_rank(1);
    if (glfwGetKey(window.window() , GLFW_KEY_2) == GLFW_PRESS) change_rank(2);
  }
}
UT_TEST_CASE_END( simplices_2d_test )

UT_TEST_CASE( simplices_3d_test )
{

}
UT_TEST_CASE_END( simplices_3d_test )

UT_TEST_CASE( simplices_4d_test )
{

}
UT_TEST_CASE_END( simplices_4d_test )

UT_TEST_CASE( polygons_test )
{

}
UT_TEST_CASE_END( polygons_test )

UT_TEST_CASE( polyhedra_test )
{

}
UT_TEST_CASE_END( polyhedra_test )


UT_TEST_SUITE_END( graphics_decomposition_suite )
