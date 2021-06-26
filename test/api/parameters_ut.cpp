//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "unit_tester.hpp"

#include "avro_params.h"

#include "adaptation/adapt.h"

#include <cmath>

using namespace avro;

UT_TEST_SUITE(api_parameters_test_suite)


UT_TEST_CASE(test1)
{

  ParameterSet parameters;

  parameters.set_param( "max parallel passes" , index_t(5) );
  UT_CATCH_EXCEPTION( parameters.set_param( "something unknown" , 2 ) );
  parameters.set_param( "metric limiting factor" , std::sqrt(2.0) );

  parameters.print();

  AdaptationParameters adapt_params(parameters);
  UT_ASSERT_EQUALS( adapt_params.get_param<index_t>("max parallel passes" ) , 5 );
  adapt_params.print();
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(api_parameters_test_suite)
