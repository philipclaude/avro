#include "unit_tester.hpp"


#include "graphics/new/tessellation.h"

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
    Field<Simplex,real_t>( topology , 3, DISCONTINUOUS )
  {
    build();

    const Simplex& element = this->element();

    for (index_t k = 0; k < topology.nb(); k++) {
      for (index_t j = 0; j < element.nb_basis(); j++) {
        // evaluate some function at the reference coordinate
        const real_t* xref = element.reference().get_reference_coordinate(j);
      }
    }
  }

};

UT_TEST_CASE( simplices_2d_test )
{
  coord_t number = 2;
  coord_t dim = number;
  std::vector<index_t> dims(number,5);
  CKF_Triangulation topology( dims );

  coord_t order = 2;
  Points nodes(dim);
  topology.element().set_basis( BasisFunctionCategory_Lagrange );
  Topology<Simplex> curvilinear(nodes,topology,order);

  std::shared_ptr<TestField> field = std::make_shared<TestField>(curvilinear);
  field->print();
  curvilinear.fields().make( "test" , field );

  Tessellation tess(number,1);
  tess.build(curvilinear);


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
