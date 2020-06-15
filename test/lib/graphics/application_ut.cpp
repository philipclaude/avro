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
#include "graphics/user_interface.h"
#include "graphics/window.h"

#include "library/ckf.h"
#include "library/obj.h"
#include "library/samples.h"

#include "mesh/topology.h"

#include "voronoi/delaunay.h"
#include "voronoi/voronoi.h"

using namespace avro;
using namespace avro::graphics;

UT_TEST_SUITE( application_suite )

UT_TEST_CASE( test1 )
{

  #if 0
  library::objFile topology( BASE_TEST_DIR+"geometry/obj/spot.obj" );
  #else
  CKF_Triangulation topology( {5,5,5} );
  #endif

  Delaunay delaunay(topology.points().dim());
  topology.points().copy(delaunay);

  Delaunay z(topology.points());
  for (index_t k=0;k<z.nb();k++)
  for (index_t d=0;d<z.dim();d++)
    z[k][d] = random_within( 0.0 , 1.0 );

  delaunay::RestrictedVoronoiDiagram rvd(topology,z);
  //rvd.compute(true);

  Visualizer vis;

  //vis.add_topology(topology);
  vis.add_topology(topology);

  // test the EPS export
  vis.run();

  vis.main_window().save_eps("tmp/test.eps");

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( application_suite )
