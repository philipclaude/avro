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

#include "shape/simplex.h"
#include "shape/quadrature.h"

using namespace avro;

UT_TEST_SUITE( shape_test_suite )

UT_TEST_CASE( simplex_shape_tests )
{
  Simplex shape( 2 , 1 );

  ConicalProductQuadrature quadrature( -1 );
  shape.load_quadrature( quadrature );

  shape.set_basis( BasisFunctionCategory_Lagrange );

  real_t x[2] = {0.5,0.5};
  real_t phi[3];
  shape.eval(x,phi);

}
UT_TEST_CASE_END( simplex_shape_tests )

UT_TEST_SUITE_END( shape_test_suite )
