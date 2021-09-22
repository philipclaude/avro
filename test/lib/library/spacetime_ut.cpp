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

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/spacetime.h"
#include "library/tesseract.h"

using namespace avro;

UT_TEST_SUITE( spacetime_topology_test_suite )

UT_TEST_CASE( test1 )
{

  std::vector<index_t> dims(4,3);
  CKF_Triangulation topology(dims);

  std::vector<real_t> x0(4,0.5);
  std::vector<real_t> length(4,1);
  library::Tesseract tesseract(x0,length);

  topology.points().attach(tesseract);

  graphics::Viewer vis;
  vis.add(topology);

  vis.run(AVRO_FULL_UNIT_TEST);
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( spacetime_topology_test_suite )
