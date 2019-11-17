#include "unit_tester.hpp"

#include "numerics/functions.h"

using namespace luna;
using namespace numerics;

UT_TEST_SUITE(FunctionsSuite)

UT_TEST_CASE(evaluation)
{
	index_t result;

	result = factorial(0);
	UT_ASSERT_EQUALS( result , 1 );

	result = factorial(1);
	UT_ASSERT_EQUALS( result , 1 );

	result = factorial(2);
	UT_ASSERT_EQUALS( result , 2 );

	result = factorial(3);
	UT_ASSERT_EQUALS( result , 6 );

	result = factorial(4);
	UT_ASSERT_EQUALS( result , 24 );

	result = factorial(5);
	UT_ASSERT_EQUALS( result , 120 );

	result = binomial(0,0);
	UT_ASSERT_EQUALS(result,1);

	result = binomial(5,2);
	UT_ASSERT_EQUALS(result,10);

	result = binomial(20,3);
	UT_ASSERT_EQUALS(result,1140);

	result = binomial(50,4);
	UT_ASSERT_EQUALS(result,230300);

	result = binomial(60,5);
	UT_ASSERT_EQUALS(result,5461512);
}
UT_TEST_CASE_END(evaluation)

UT_TEST_SUITE_END(FunctionsSuite)
