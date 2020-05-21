#include "unit_tester.hpp"

#include "master/simplex.h"
#include "master/quadrature.h"

using namespace avro;

UT_TEST_SUITE( master_test_suite )

UT_TEST_CASE( simplex_master_tests )
{
  Simplex master( 2 , 1 );

  ConicalProductQuadrature quadrature( -1 );
  master.load_quadrature( quadrature );

  master.set_basis( BasisFunctionCategory_Lagrange );

  real_t x[2] = {0.5,0.5};
  real_t phi[3];
  master.eval(x,phi);

}
UT_TEST_CASE_END( simplex_master_tests )

UT_TEST_SUITE_END( master_test_suite )
