// luna: unstrucutred adaptation library
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "unit_tester.hpp"

#include "common/types.h"

#include "numerics/vector.h"

using namespace luna;
using namespace luna::numerics;

UT_TEST_SUITE(VectorSuite)

UT_TEST_CASE(test1)
{
  Vector<real_t> x( {1,2,3} );
  Vector<real_t> y( {4,5,6} );

  x = x +y;

  UT_ASSERT_EQUALS( x[0] , 5 );
  UT_ASSERT_EQUALS( x[1] , 7 );
  UT_ASSERT_EQUALS( x[2] , 9 );

  Vector<real_t> z = y*3;
  UT_ASSERT_EQUALS( z[0] , 12 );
  UT_ASSERT_EQUALS( z[1] , 15 );
  UT_ASSERT_EQUALS( z[2] , 18 );


}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(VectorSuite)
