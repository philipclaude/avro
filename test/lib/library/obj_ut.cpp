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

#include "graphics/application.h"
#include "library/obj.h"

using namespace avro;
using namespace avro::library;

UT_TEST_SUITE( obj_file_suite )

UT_TEST_CASE( test1 )
{

  objFile topology( BASE_TEST_DIR+"/geometry/obj/suzanne.obj" );
  //objFile topology( "/Users/pcaplan/Google Drive/teaching/cs461w21/assignments/project/giraffe/model/giraffe.obj" );

  graphics::Visualizer vis;
  vis.add_topology(topology);
  vis.run();

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( obj_file_suite )
