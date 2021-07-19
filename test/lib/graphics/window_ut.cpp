#include "unit_tester.hpp"

#include "graphics/managers.h"
#include "graphics/vao.h"
#include "graphics/window.h"

#include "graphics/colormap.h"
#include "graphics/shader.h"

#include "library/ckf.h"
#include "library/egads.h"

#include "mesh/field.hpp"
#include "mesh/points.h"
#include "mesh/topology.h"

#if AVRO_HEADLESS_GRAPHICS == 0

using namespace avro;
using namespace avro::graphics;

UT_TEST_SUITE( graphics_window_suite )

class TestField : public Field<Simplex,std::vector<real_t>> {

public:
  TestField( const Topology<Simplex>& topology , coord_t order ) :
    Field( topology , order, DISCONTINUOUS )
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

UT_TEST_CASE( simplices_2d_test )
{
  coord_t number = 3;
  coord_t dim = number;

  EGADS::Context context;
  std::vector<real_t> lengths(number,1.0);
  EGADS::Cube geometry(&context,lengths);

  std::vector<index_t> dims(number,3);
  //dims[1] = 2;
  CKF_Triangulation topology( dims );
  topology.element().set_basis( BasisFunctionCategory_Lagrange );
  topology.points().attach(geometry);

  coord_t geometry_order = 2;
  Points nodes(dim);
  topology.element().set_basis( BasisFunctionCategory_Lagrange );
  Topology<Simplex> curvilinear(nodes,topology,geometry_order);
  curvilinear.element().set_basis( BasisFunctionCategory_Lagrange );
  curvilinear.points().attach(geometry);

  #if 1
  for (index_t k = 0; k < curvilinear.points().nb(); k++) {

    real_t s = curvilinear.points()[k][0];
    real_t t = curvilinear.points()[k][1];

    real_t theta = s*M_PI;
    real_t R = 0.75 + t*(1.0 - 0.75);

    curvilinear.points()[k][0] = R*cos(M_PI - theta);
    curvilinear.points()[k][1] = R*sin(theta);
  }
  #endif

  if (AVRO_FULL_UNIT_TEST) return;

  coord_t solution_order = 3;
  std::shared_ptr<TestField> field = std::make_shared<TestField>(curvilinear,solution_order);
  field->element().set_basis( BasisFunctionCategory_Lagrange );
  curvilinear.fields().make( "test" , field );

  int width = 400, height = width;
  Window window(width,height);
  window.init();

  Plot plot(curvilinear);
  plot.build();
  window.add_plot(&plot);

  Plot plot2(topology);
  plot2.build();
  window.add_plot(&plot2);

  window.compute_view();
  window.enable_controls(true);

  // bind the colormap values to a buffer
  gl_index colormap_buffer;
  GL_CALL( glGenBuffers( 1 , &colormap_buffer ) );
  Colormap colormap;
  colormap.change_style("giraffe");
  index_t ncolor = 256*3;
  GL_CALL( glBindBuffer( GL_TEXTURE_BUFFER , colormap_buffer) );
  GL_CALL( glBufferData( GL_TEXTURE_BUFFER , sizeof(float) * ncolor , colormap.data() , GL_STATIC_DRAW) );

  // generate a texture to hold the colormap buffer
  GLuint colormap_texture;
  GL_CALL( glGenTextures( 1 , &colormap_texture) );
  GL_CALL( glActiveTexture( GL_TEXTURE0 + 1 ) );
  GL_CALL( glBindTexture( GL_TEXTURE_BUFFER , colormap_texture) );
  GL_CALL( glTexBuffer( GL_TEXTURE_BUFFER , GL_R32F , colormap_buffer ) );

  // initial draw, subsequent drawing will only be performed when a callback is invoked
  window.draw();
  while (true) {

    window.draw();

    // wait for user input
    glfwWaitEvents();

    // determine if we should exit the render loop
    if (glfwWindowShouldClose(window.window())) break;
    if (glfwGetKey(window.window(), GLFW_KEY_ESCAPE ) == GLFW_PRESS) break;
  }
}
UT_TEST_CASE_END( simplices_2d_test )

UT_TEST_SUITE_END( graphics_window_suite )

#else // AVRO_HEADLESS_GRAPHICS is equal to 1
UT_TEST_SUITE( graphics_decomposition_suite )
UT_TEST_SUITE_END( graphics_decomposition_suite )
#endif
