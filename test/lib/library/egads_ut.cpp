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

#include "geometry/egads/context.h"

#include "library/egads.h"

using namespace avro;

UT_TEST_SUITE( library_egads_suite )

UT_TEST_CASE( box_test )
{
  EGADS::Context context;
  real_t x0[3] = {0,0,0};
  EGADS::Cube cube(&context,{1,1,1},x0);
}
UT_TEST_CASE_END( box_test )

UT_TEST_SUITE_END( library_egads_suite )
