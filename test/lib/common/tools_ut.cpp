#include "unit_tester.hpp"

#include "common/tools.h"
#include "common/types.h"

#include <string>
#include <vector>

using namespace avro;

UT_TEST_SUITE(tools_test_suite)


UT_TEST_CASE(test1)
{
  std::vector<index_t> f = {1,2,3};
  std::string label = unique_label<index_t>(f);
  UT_ASSERT_EQUALS( label , "1|2|3" );

  std::string txt = "3.2";
  real_t val_float = unstringify<real_t>(txt);
  UT_ASSERT_NEAR( val_float , 3.2 , 1e-12 );

  txt = "6";
  UT_ASSERT_EQUALS( unstringify<index_t>(txt) , 6 );

  txt = "True";
  UT_ASSERT_EQUALS( unstringify<bool>(txt) , true );


  index_t x = 3;
  printValue(x);

  printValue<float>(3.14);
  printValue<double>(3.14);

  real_t r0 = random_within( 0.0 , 1.0 );
  UT_ASSERT( r0 >= 0.0 and r0 <= 1.0 );

  int r1 = random_within( -100 , 100 );
  UT_ASSERT( r1 >= -100 and r1 <= 100 );
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(tools_test_suite)
