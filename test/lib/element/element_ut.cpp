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

#include "element/simplex.h"
#include "element/quadrature.h"

using namespace avro;

UT_TEST_SUITE( element_test_suite )

UT_TEST_CASE( simplex_element_tests )
{
  Simplex element( 2 , 1 );

  ConicalProductQuadrature quadrature( -1 );
  element.load_quadrature( quadrature );

  element.set_basis( BasisFunctionCategory_Lagrange );

  real_t x[2] = {0.5,0.5};
  real_t phi[3];
  element.eval(x,phi);

}
UT_TEST_CASE_END( simplex_element_tests )

UT_TEST_SUITE_END( element_test_suite )
