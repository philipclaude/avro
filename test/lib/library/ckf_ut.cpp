#include "unit_tester.hpp"

#include "library/ckf.h"

using namespace luna;

UT_TEST_SUITE( ckf_suite )

UT_TEST_CASE( test_2d )
{

  index_t n = 50;

  CKF_Triangulation topology( {n,n,n} );

  printf("generated %lu %lu-simplices\n",topology.nb(),topology.number());

}
UT_TEST_CASE_END( test_2d )

UT_TEST_SUITE_END( ckf_suite )
