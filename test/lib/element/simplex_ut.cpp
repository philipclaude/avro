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

#include "element/simplex.h"
#include "element/quadrature.h"

using namespace avro;

UT_TEST_SUITE( simplex_test_suite )

UT_TEST_CASE( simplex_element_tests )
{
  coord_t number = 1;
  coord_t order  = 1;

  Simplex element( number, order );

}
UT_TEST_CASE_END( simplex_element_tests )

UT_TEST_SUITE_END( simplex_test_suite )
