//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
// avro: unstrucutred adaptation library
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "unit_tester.hpp"

#include "avro_types.h"

#include "numerics/linear_algebra.h"
#include "numerics/vec.h"

using namespace avro;
using namespace avro::numerics;

UT_TEST_SUITE(vector_suite)

UT_TEST_CASE(test1)
{
  vecd<real_t> x( {1,2,3} );
  vecd<real_t> y( {4,5,6} );

  x = x +y;

  UT_ASSERT_EQUALS( x[0] , 5 );
  UT_ASSERT_EQUALS( x[1] , 7 );
  UT_ASSERT_EQUALS( x[2] , 9 );

  /*vecd<real_t> z = 3.0*y;
  UT_ASSERT_EQUALS( z[0] , 12 );
  UT_ASSERT_EQUALS( z[1] , 15 );
  UT_ASSERT_EQUALS( z[2] , 18 );*/


}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(vector_suite)
