#include "unit_tester.hpp"

#include "graphics/new/managers.h"
#include "graphics/new/vao.h"
#include "graphics/new/application.h"

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

UT_TEST_SUITE( graphics_decomposition_suite )

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

  std::vector<index_t> dims(number,10);
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

  coord_t solution_order = 3;
  std::shared_ptr<TestField> field = std::make_shared<TestField>(curvilinear,solution_order);
  field->element().set_basis( BasisFunctionCategory_Lagrange );
  curvilinear.fields().make( "test" , field );

  OpenGL_Application app;
  app.add( curvilinear );
  app.add( topology );

  app.run();

}
UT_TEST_CASE_END( simplices_2d_test )

UT_TEST_SUITE_END( graphics_decomposition_suite )

#else // AVRO_HEADLESS_GRAPHICS is equal to 1
UT_TEST_SUITE( graphics_decomposition_suite )
UT_TEST_SUITE_END( graphics_decomposition_suite )
#endif
