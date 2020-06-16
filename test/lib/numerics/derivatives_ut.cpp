#include "unit_tester.hpp"

#include "numerics/dual.h"

#include <numpack/types/SurrealS.h>

using namespace avro;

template<typename type>
type
test_function( type& x )
{
  return x*sin(x) + exp(x);
}

real_t
test_function_derivative( real_t x )
{
  return sin(x) + x*cos(x) + exp(x);
}

UT_TEST_SUITE( derivatives_test_suite )

UT_TEST_CASE( dual_tests )
{
  real_t x = 1.2;
  dual xd(x,1); // set the dual component of the variable to 1 to comppute the derivative

  std::cout << xd << std::endl;

  dual fd = test_function(xd);
  real_t df_dx = fd.du();
  real_t df_dx0 = test_function_derivative(x);

  printf("derivative = %g\n",df_dx);
  printf("analytic = %g\n",df_dx0);

  UT_ASSERT_NEAR( df_dx , df_dx0 , 1e-12 );
}
UT_TEST_CASE_END( dual_tests )

UT_TEST_CASE( surreal_tests )
{

  real_t x = 1.2;

  SurrealS<1> xs(x);
  xs.deriv(0) = 1.0;

  SurrealS<1> fs = test_function(xs);
  std::cout << fs << std::endl;

  real_t df_dx = fs.deriv(0);
  real_t df_dx0 = test_function_derivative(x);

  UT_ASSERT_NEAR( df_dx , df_dx0 , 1e-12 );

}
UT_TEST_CASE_END( surreal_tests )

UT_TEST_SUITE_END( derivatives_test_suite )
