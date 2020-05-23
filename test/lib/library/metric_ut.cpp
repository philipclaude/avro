#include "unit_tester.hpp"

#include "library/metric.h"

using namespace avro;
using namespace avro::library;

UT_TEST_SUITE( library_metric_test_suite )

UT_TEST_CASE( test1 )
{
  // not much to test here i guess...

  MetricField_Uniform uniform1(3,{1,2,3});
  MetricField_Uniform uniform2(3,1);
  MetricField_UGAWG_Linear linear;
  MetricField_UGAWG_Polar1 polar1;
  MetricField_UGAWG_Polar2 polar2;
  MetricField_Tesseract_Linear tl;
  MetricField_Tesseract_Wave tw;

  numerics::SymMatrixD<real_t> m;
  real_t x3[3] = {1,1,1};

  m = uniform1(x3);
  m = uniform2(x3);
  m = linear(x3);
  m = polar1(x3);
  m = polar2(x3);

  real_t x4[4] {1,1,1,1};
  m = tl(x4);
  m = tw(x4);
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( library_metric_test_suite )
