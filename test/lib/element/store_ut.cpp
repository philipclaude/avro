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

#include "element/quadrature.h"

using namespace avro;

UT_TEST_SUITE( quadrature_cache_test_suite )

UT_TEST_CASE( build )
{
  const matd<real_t>& B = __store_simplex_lagrange__.get_basis( 3 , 3 , 4 );

  B.print();

  const Quadrature& quadrature = __store_simplex_lagrange__.quadrature( 3 , 4 );

  printf("[status] memory usage for quadrature store: %lu kB\n",__store_simplex_lagrange__.memory()/1000);

}
UT_TEST_CASE_END( build )

UT_TEST_SUITE_END( quadrature_cache_test_suite )
