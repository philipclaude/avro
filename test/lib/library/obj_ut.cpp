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
#include "library/obj.h"

using namespace avro;
using namespace avro::library;

UT_TEST_SUITE( obj_file_suite )

UT_TEST_CASE( test1 )
{

  //objFile topology( AVRO_SOURCE_DIR+"/build/release/bunny.obj" );

  objFile topology( BASE_TEST_DIR+"/geometry/obj/suzanne.obj" );
  //objFile topology( "/Users/pcaplan/Google Drive/teaching/cs461w21/assignments/project/giraffe/model/giraffe.obj" );

  graphics::Viewer vis;
  vis.add(topology);
  vis.run(AVRO_FULL_UNIT_TEST);

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( obj_file_suite )
