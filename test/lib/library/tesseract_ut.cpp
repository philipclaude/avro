#include "unit_tester.hpp"

#include "library/tesseract.h"

using namespace avro;

UT_TEST_SUITE( tesseract_suite )

UT_TEST_CASE( test1 )
{
  std::vector<real_t> x0(4,0);
  std::vector<real_t> length(4,1);

  library::Tesseract tesseract(x0,length);

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( tesseract_suite )
