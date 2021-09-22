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

#include "avro_params.h"

#include "adaptation/adapt.h"

#include <cmath>

using namespace avro;

UT_TEST_SUITE(api_parameters_test_suite)


UT_TEST_CASE(test1)
{

  ParameterSet parameters;

  parameters.set( "max parallel passes" , index_t(5) );
  UT_CATCH_EXCEPTION( parameters.set( "something unknown" , index_t(2)  ) );
  parameters.set( "metric limiting factor" , std::sqrt(2.0) );
  parameters.set( "adapt iter" , index_t(10) );

  parameters.print();

  AdaptationParameters adapt_params(parameters);
  index_t max_parallel_passes = adapt_params["max parallel passes"];
  UT_ASSERT_EQUALS( max_parallel_passes  , 5 );
  adapt_params.print();

  // even though the default is 1, the inherited parameter should not be overwritten
  index_t adapt_iter = adapt_params["adapt iter"];
  UT_ASSERT_EQUALS( adapt_iter , 10 );
  int adapt_iter_int;
  UT_CATCH_EXCEPTION( adapt_iter_int = adapt_params["adapt iter"] );

  index_t ivf;
  UT_CATCH_EXCEPTION( ivf = adapt_params["insertion volume factor"] );

  std::string directory = adapt_params["directory"];
  printf("dir = %s\n",directory.c_str());

  bool swapout = adapt_params["swapout"];
  UT_ASSERT_EQUALS( swapout , false );
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(api_parameters_test_suite)
