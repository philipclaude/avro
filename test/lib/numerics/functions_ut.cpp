//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "unit_tester.hpp"

#include "numerics/functions.h"

using namespace avro;
using namespace numerics;

UT_TEST_SUITE(functions_suite)

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

UT_TEST_SUITE_END(functions_suite)
