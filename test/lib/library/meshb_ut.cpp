#include "unit_tester.hpp"

#include "geometry/egads/context.h"
#include "geometry/egads/model.h"

#include "library/meshb.h"

using namespace luna;

UT_TEST_SUITE( meshb_test_suite )

UT_TEST_CASE( test1 )
{
  EGADS::Context context;
  EGADS::Model model(context,"/Users/pcaplan/Desktop/cube-cylinder.egads");

  library::meshb mesh( "/Users/pcaplan/Desktop/cube-cylinder.mesh" , &model );

  mesh.points().compute_params();

  mesh.points().print(true);
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( meshb_test_suite )
