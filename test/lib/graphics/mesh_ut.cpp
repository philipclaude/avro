#include "unit_tester.hpp"

#include "graphics/new/tessellation.h"
#include "graphics/new/window.h"

#include "graphics/shader.h"

#include "library/ckf.h"

#include "mesh/field.hpp"
#include "mesh/points.h"
#include "mesh/topology.h"

using namespace avro;
using namespace avro::graphics;

UT_TEST_SUITE( graphics_decomposition_suite )

class TestField : public Field<Simplex,real_t> {

public:
  TestField( const Topology<Simplex>& topology ) :
    Field<Simplex,real_t>( topology , 2, DISCONTINUOUS )
  {
    build();

    const Simplex& element = this->element();

    for (index_t k = 0; k < topology.nb(); k++) {
      for (index_t j = 0; j < element.nb_basis(); j++) {
        // evaluate some function at the reference coordinate
        const real_t* xref = element.reference().get_reference_coordinate(j);
        UNUSED(xref);
      }
    }
  }
};

graphics::mat4 model_matrix = glm::identity();

index_t width = 400;
index_t height = 400;
bool dragging = false;
int xm = -1;
int ym = -1;
int XM,YM;
Tessellation* tess_ptr = nullptr;
ShaderProgram* triangle_shader = nullptr;
ShaderProgram* edge_shader = nullptr;

static void
mouse_button_callback(GLFWwindow* window ,int button,int action,int mods)
{
  if (action == GLFW_PRESS) dragging = true;
  else {
    xm = ym = -1;
    dragging = false;
  }
}

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

void
draw() {
  float col = 1.0;
  glClearColor (col,col,col, 0.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  tess_ptr->draw_triangles(*triangle_shader);
  tess_ptr->draw_edges(*edge_shader);
}

static void
mouse_move_callback(GLFWwindow* window, double x, double y)
{
  if (!dragging) return;
  if (xm < 0 || ym < 0) {
    xm = x;
    ym = y;
    return;
  }

  // for now, hack the center of rotation to be the center of the cube domain
  real_t xc = 0.5;
  real_t yc = 0.5;
  real_t zc = 0.5;
  mat4 T = glm::translate( glm::identity() , {xc,yc,zc} );
  mat4 Tinv = glm::translate( glm::identity() , {-xc,-yc,-zc} );

  // translate to the center, rotate, then translate back
  mat4 R = rotation( -(x - xm)/width , (y - ym)/height );
  model_matrix = T * R * Tinv * model_matrix;

  xm = x;
  ym = y;

  draw();
}

UT_TEST_CASE( simplices_2d_test )
{
  coord_t number = 3;
  coord_t dim = number;
  std::vector<index_t> dims(number,3);
  CKF_Triangulation topology( dims );

  coord_t order = 2;
  Points nodes(dim);
  topology.element().set_basis( BasisFunctionCategory_Lagrange );
  Topology<Simplex> curvilinear(nodes,topology,order);

  for (index_t k = 0; k < curvilinear.points().nb(); k++) {

    real_t s = curvilinear.points()[k][0];
    real_t t = curvilinear.points()[k][1];

    real_t theta = s*M_PI;
    real_t R = 0.5 + t*(1.0 - 0.5);

    curvilinear.points()[k][0] = R*cos(theta);
    curvilinear.points()[k][1] = R*sin(theta);
  }

  std::shared_ptr<TestField> field = std::make_shared<TestField>(curvilinear);
  field->print();
  curvilinear.fields().make( "test" , field );

  Window window(width,height);
  window.init();

  Tessellation tess(number,1);
  tess.build(curvilinear);

  bool with_tess = true;
  ShaderProgram tshader("triangles",with_tess);
  tshader.use();

  ShaderProgram eshader("edges",with_tess);
  eshader.use();

  // set up the view
  graphics::vec3 eye = {0,0,5};
  graphics::vec3 up = {0,1,0};
  graphics::vec3 center = {0.5,0.5,0.0};

  float fov    = M_PI/4.0;
  float aspect = width/height;
  float znear  = 0.001;
  float zfar   = 100;

  // define the transformation matrices, and set the MVP into the shader program
  graphics::mat4 perspective_matrix = glm::perspective( fov , aspect , znear , zfar );
  graphics::mat4 view_matrix = glm::lookAt( eye , center , up );
  graphics::mat4 mv = view_matrix * model_matrix;
  graphics::mat4 mvp = perspective_matrix * mv;

  tshader.use();
  tshader.setUniform("u_ModelViewProjectionMatrix",mvp);

  eshader.use();
  eshader.setUniform("u_ModelViewProjectionMatrix",mvp);

  glfwSetCursorPosCallback(window.window(),&mouse_move_callback);
  glfwSetMouseButtonCallback(window.window(),&mouse_button_callback);

  // set the tessellation level for the TCS
  int level = 10;

  tshader.use();
  tshader.setUniform( "nb_basis" , 6 );
  tshader.setUniform( "u_level" , level );

  eshader.use();
  eshader.setUniform("nb_basis" , 3 );
  eshader.setUniform("u_level" , level );

  tess_ptr = &tess;
  triangle_shader = &tshader;
  edge_shader = &eshader;
  draw();

  while (true) {

    // todo, only set the shader uniforms when drawing
    graphics::mat4 mv = view_matrix * model_matrix;
    graphics::mat4 mvp = perspective_matrix * mv;

    tshader.use();
    tshader.setUniform("u_ModelViewProjectionMatrix",mvp);

    eshader.use();
    eshader.setUniform("u_ModelViewProjectionMatrix",mvp);

    // swap buffers and wait for user input
    glfwSwapBuffers(window.window());
    glfwWaitEvents();
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
