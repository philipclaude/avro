#include "unit_tester.hpp"

#include "numerics/dual.h"

using namespace avro;

UT_TEST_SUITE(dual_test_suite)

UT_TEST_CASE(operations)
{
  double r=2.2;
  double tol=1e-12;
  dual u(2.5,3.6);
  dual v(1.89,4.02);
  dual w,z(r);
  bool logical;

  /* Assignment and initialization */
  UT_ASSERT_NEAR(2.500000000000000e+00,u.re(),tol);
  UT_ASSERT_NEAR(3.600000000000000e+00,u.du(),tol);
  UT_ASSERT_NEAR(0.000000000000000e+00,w.re(),tol);
  UT_ASSERT_NEAR(0.000000000000000e+00,w.du(),tol);
  UT_ASSERT_NEAR(2.200000000000000e+00,z.re(),tol);
  UT_ASSERT_NEAR(0.000000000000000e+00,z.du(),tol);
  w = u;
  UT_ASSERT_NEAR(2.500000000000000e+00,w.re(),tol);
  UT_ASSERT_NEAR(3.600000000000000e+00,w.du(),tol);

  /* Addition */
  w = u +v;
  UT_ASSERT_NEAR(4.390000000000000e+00,w.re(),tol);
  UT_ASSERT_NEAR(7.619999999999999e+00,w.du(),tol);
  w = u +r;
  UT_ASSERT_NEAR(4.700000000000000e+00,w.re(),tol);
  UT_ASSERT_NEAR(3.600000000000000e+00,w.du(),tol);
  w = r +u;
  UT_ASSERT_NEAR(4.700000000000000e+00,w.re(),tol);
  UT_ASSERT_NEAR(3.600000000000000e+00,w.du(),tol);

  /* Subtraction */
  w = u -v;
  UT_ASSERT_NEAR(6.100000000000001e-01,w.re(),tol);
  UT_ASSERT_NEAR(-4.199999999999995e-01,w.du(),tol);
  w = u -r;
  UT_ASSERT_NEAR(2.999999999999998e-01,w.re(),tol);
  UT_ASSERT_NEAR(3.600000000000000e+00,w.du(),tol);
  w = r -u;
  UT_ASSERT_NEAR(-2.999999999999998e-01,w.re(),tol);
  UT_ASSERT_NEAR(-3.600000000000000e+00,w.du(),tol);

  /* Multiplication */
  w = u*v;
  UT_ASSERT_NEAR(4.725000000000000e+00,w.re(),tol);
  UT_ASSERT_NEAR(1.685400000000000e+01,w.du(),tol);
  w = r*u;
  UT_ASSERT_NEAR(5.500000000000000e+00,w.re(),tol);
  UT_ASSERT_NEAR(7.920000000000001e+00,w.du(),tol);
  w = u*r;
  UT_ASSERT_NEAR(5.500000000000000e+00,w.re(),tol);
  UT_ASSERT_NEAR(7.920000000000001e+00,w.du(),tol);

  /* Division */
  w = u/v;
  UT_ASSERT_NEAR(1.322751322751323e+00,w.re(),tol);
  UT_ASSERT_NEAR(-9.087091626774162e-01,w.du(),tol);
  w = r/u;
  UT_ASSERT_NEAR(8.800000000000001e-01,w.re(),tol);
  UT_ASSERT_NEAR(-1.267200000000000e+00,w.du(),tol);
  w = u/r;
  UT_ASSERT_NEAR(1.136363636363636e+00,w.re(),tol);
  UT_ASSERT_NEAR(1.636363636363636e+00,w.du(),tol);

  /* Math */
  w = pow(u,r);
  UT_ASSERT_NEAR(7.507027712383946e+00,w.re(),tol);
  UT_ASSERT_NEAR(2.378226379283235e+01,w.du(),tol);
  w = pow(u,v);
  UT_ASSERT_NEAR(5.650756800854586e+00,w.re(),tol);
  UT_ASSERT_NEAR(3.619359876969143e+01,w.du(),tol);
  w = exp(u);
  UT_ASSERT_NEAR(1.218249396070347e+01,w.re(),tol);
  UT_ASSERT_NEAR(4.385697825853251e+01,w.du(),tol);
  w = log(u);
  UT_ASSERT_NEAR(9.162907318741551e-01,w.re(),tol);
  UT_ASSERT_NEAR(1.440000000000000e+00,w.du(),tol);
  w = sqrt(u);
  UT_ASSERT_NEAR(1.581138830084190e+00,w.re(),tol);
  UT_ASSERT_NEAR(1.138419957660616e+00,w.du(),tol);
  w = sin(u);
  UT_ASSERT_NEAR(5.984721441039565e-01,w.re(),tol);
  UT_ASSERT_NEAR(-2.884117015968962e+00,w.du(),tol);
  w = cos(u);
  UT_ASSERT_NEAR(-8.011436155469337e-01,w.re(),tol);
  UT_ASSERT_NEAR(-2.154499718774244e+00,w.du(),tol);
  w = tan(u);
  UT_ASSERT_NEAR(-7.470222972386603e-01,w.re(),tol);
  UT_ASSERT_NEAR(5.608952325258212e+00,w.du(),tol);

  /* Logicals */
  logical = u >v;
  UT_ASSERT_EQUALS(true,logical);
  logical = r >u;
  UT_ASSERT_EQUALS(false,logical);
  logical = u >r;
  UT_ASSERT_EQUALS(true,logical);
  logical = u <v;
  UT_ASSERT_EQUALS(false,logical);
  logical = r <u;
  UT_ASSERT_EQUALS(true,logical);
  logical = u <r;
  UT_ASSERT_EQUALS(false,logical);

  /* Printing */
  u.print();
  std::cout << v << std::endl;

  return;
}
UT_TEST_CASE_END(operations)

UT_TEST_SUITE_END(dual_test_suite)
