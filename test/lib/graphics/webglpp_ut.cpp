#include "unit_tester.hpp"

#include "graphics/application.h"
#include "graphics/vao.h"
#include "graphics/window.h"

#include "graphics/colormap.h"
#include "graphics/shader.h"
#include "graphics/webglpp.h"

#include "library/ckf.h"
#include "library/egads.h"

#include "mesh/field.hpp"
#include "mesh/points.h"
#include "mesh/topology.h"

using namespace avro;
using namespace avro::graphics;

class TestField : public Field<Simplex,real_t> {

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

        x[3] = 0.95*sin( 10*x[0]*M_PI )*sin( 10*x[1]*M_PI ) +
               0.95*sin( 10*x[1]*M_PI )*sin( 10*x[2]*M_PI );
        if (x[3] < 0) x[3] = 0;
        this->value(idx) = x[3];
      }
    }
  }

  index_t nb_rank() const { return 1; }
};

UT_TEST_SUITE( webglpp_test_suite )

UT_TEST_CASE(test1)
{
  coord_t number = 3;
  coord_t dim = number;
  std::vector<index_t> dims(number,50);
  CKF_Triangulation topology( dims );

  coord_t geometry_order = 1;
  Points nodes(dim);
  topology.element().set_basis( BasisFunctionCategory_Lagrange );
  Topology<Simplex> curvilinear(nodes,topology,geometry_order);
  curvilinear.element().set_basis( BasisFunctionCategory_Lagrange );

  EGADS::Context context;
  std::vector<real_t> lengths(number,1.0);
  EGADS::Cube geometry(&context,lengths);
  curvilinear.points().attach(geometry);

  coord_t solution_order = 3;
  std::shared_ptr<TestField> field = std::make_shared<TestField>(curvilinear,solution_order);
  field->element().set_basis( BasisFunctionCategory_Lagrange );
  curvilinear.fields().make( "test" , field );


  Viewer viewer(false);
  viewer.add(curvilinear);
  viewer.run();
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END( webglpp_test_suite )
