#include "unit_tester.hpp"

#include "geometry/egads/context.h"

#include "library/egads.h"

using namespace luna;

UT_TEST_SUITE( library_egads_suite )

UT_TEST_CASE( box_test )
{
  EGADS::Context context;
  real_t x0[3] = {0,0,0};
  EGADS::Cube cube(&context,{1,1,1},x0);
}
UT_TEST_CASE_END( box_test )

UT_TEST_SUITE_END( library_egads_suite )
