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
  curvilinear.Table<index_t>::print();

  std::shared_ptr<TestField> field = std::make_shared<TestField>(curvilinear);
  field->print();
  curvilinear.fields().make( "test" , field );

  index_t width = 400;
  index_t height = 400;
  Window window(width,height);
  window.init();

  Tessellation tess(number,1);
  tess.build(curvilinear);

  bool with_tess = true;
  ShaderProgram shader("basic",with_tess);
  shader.use();

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
  graphics::mat4 model_matrix = glm::identity();
  graphics::mat4 mv = view_matrix * model_matrix;
  graphics::mat4 mvp = perspective_matrix * mv;
  shader.setUniform("u_ModelViewProjectionMatrix",mvp);

  int level = 1;
  shader.setUniform( "nb_basis" , 10 );
  shader.setUniform( "u_level" , level );

  while (true) {

    // grey-ish background color
    float col = 0.9;
    glClearColor (col,col,col, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    tess.draw();

    // swap buffers and wait for user input
    glfwSwapBuffers(window.window());
    glfwPollEvents();

    //break;
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
