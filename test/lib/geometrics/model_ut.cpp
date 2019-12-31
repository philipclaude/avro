#include "unit_tester.hpp"

#include "geometry/egads/context.h"
#include "geometry/egads/model.h"

using namespace avro;

UT_TEST_SUITE( geometry_egads_model )

UT_TEST_CASE( test_io )
{
  EGADS::Context context;
  EGADS::Model model(context,"/Users/pcaplan/Desktop/cube-cylinder.egads" );

  model.body(0).print();

}
UT_TEST_CASE_END( test_io )

UT_TEST_SUITE_END( geometry_egads_model )
