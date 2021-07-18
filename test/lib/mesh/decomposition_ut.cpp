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

#include "library/ckf.h"
#include "library/samples.h"

#include "mesh/decomposition.h"

#include "voronoi/delaunay.h"
#include "voronoi/voronoi.h"

using namespace avro;
using namespace avro::graphics;

UT_TEST_SUITE( triangulation_suite )

UT_TEST_CASE( simplex_topologies )
{

}
UT_TEST_CASE_END( simplex_topologies )

UT_TEST_CASE( voronoi_tests )
{
  CKF_Triangulation topology( {4,4,4} );
  library::RegularPolygon polygon(6);

  Delaunay z(topology.points());
  for (index_t k=0;k<z.nb();k++)
  for (index_t d=0;d<z.dim();d++)
    z[k][d] = random_within( 0.0 , 1.0 );

  delaunay::RestrictedVoronoiDiagram rvd(topology,z);
  rvd.compute();

  graphics::Viewer vis;
  vis.add(rvd);
  //vis.run();
}
UT_TEST_CASE_END( voronoi_tests )

UT_TEST_SUITE_END( triangulation_suite )
