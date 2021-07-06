#include "unit_tester.hpp"

#include "element/quadrature.h"

using namespace avro;

UT_TEST_SUITE( quadrature_cache_test_suite )

UT_TEST_CASE( build )
{
  const matd<real_t>& B = __store_simplex_lagrange__.get_basis( 3 , 3 , 5 );

  B.print();

  const Quadrature& quadrature = __store_simplex_lagrange__.quadrature( 3 , 7 );

  printf("[status] memory usage for quadrature store: %lu kB\n",__store_simplex_lagrange__.memory()/1000);

}
UT_TEST_CASE_END( build )

UT_TEST_SUITE_END( quadrature_cache_test_suite )
