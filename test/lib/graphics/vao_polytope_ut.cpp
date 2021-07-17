#include "unit_tester.hpp"

#include "graphics/new/managers.h"
#include "graphics/new/vao.h"
#include "graphics/new/window.h"

#include "graphics/colormap.h"
#include "graphics/shader.h"

#include "library/ckf.h"
#include "library/egads.h"

#include "mesh/points.h"
#include "mesh/topology.h"

#include "voronoi/delaunay.h"
#include "voronoi/optimal_transport.h"

#if AVRO_HEADLESS_GRAPHICS == 0

using namespace avro;
using namespace avro::graphics;

UT_TEST_SUITE( graphics_decomposition_suite )

class TestField : public Field<Polytope,real_t> {

public:
  TestField( Topology<Polytope>& topology ) :
    Field( topology , 0, DISCONTINUOUS )
  {
    build();

    index_t n = 0;
    for (index_t k = 0; k < topology.nb(); k++) {
      this->value(n++) = k;
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
graphics::mat4 normal_matrix = glm::transpose(glm::inverse(mv));

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
  normal_matrix = glm::transpose(glm::inverse(mv));

  triangle_shader->use();
  triangle_shader->setUniform("u_ModelViewProjectionMatrix",mvp);
  triangle_shader->setUniform("u_NormalMatrix",normal_matrix);
  triangle_shader->setUniform("u_ModelViewMatrix",mv);

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
mouse_button_callback(GLFWwindow* window , int button, int action, int mods) {
  if (action == GLFW_PRESS) dragging = true;
  else {
    xm = ym = -1;
    dragging = false;
  }
}

static void
mouse_move_callback(GLFWwindow* window, double x, double y) {
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

UT_TEST_CASE( vao_polytopes_test )
{
  typedef Polytope type;
  coord_t number = 2;
  coord_t dim = number;
  index_t nb_points = 1e2;

  // create random delaunay vertices
  CubeDomain<type> domain(dim,dim,10);
  Delaunay delaunay( dim );
  std::vector<index_t> elems;
  for (index_t k = 0; k < nb_points; k++) {

    index_t elem = 0;
    std::vector<real_t> p(dim,0.);
    for (coord_t d = 0; d < dim; d++)
      p[d] = random_within(0.0,1.0);

    // create the delaunay point and retain which element in the domain it is in
    delaunay.create(p.data());
    elems.push_back( elem );
  }

  // initialize and compute the laguerre diagram
  delaunay::LaguerreDiagram<type> diagram( delaunay , domain );
  diagram.set_elements( elems );

  delaunay::IntegrationSimplices triangulation(number,number);
  diagram.compute(false,&triangulation);

  std::shared_ptr<TestField> fld = std::make_shared<TestField>(diagram);
  diagram.fields().make("test",fld);


  Window window(width,height);
  window.init();
  win = window.window();

  VertexAttributeObject vao;
  vao.build(diagram);
  window.manager().write(vao);

  std::vector<std::string> macros = {"#define WITH_TESSELLATION 0","#define SOLUTION_ORDER -1",
                                     "#define GEOMETRY_ORDER 1"};

  bool with_tess = false;
  ShaderProgram tshader("triangles",with_tess,macros);
  ShaderProgram eshader("edges",with_tess,macros);
  ShaderProgram pshader("points");

  tshader.setUniform("u_umin",0);
  tshader.setUniform("u_umax",300);

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
  int level = 1;

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
  }
}
UT_TEST_CASE_END( vao_polytopes_test )


UT_TEST_SUITE_END( graphics_decomposition_suite )

#else
UT_TEST_SUITE( graphics_decomposition_suite )
UT_TEST_SUITE_END( graphics_decomposition_suite )

#endif
