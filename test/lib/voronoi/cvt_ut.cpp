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

#include "voronoi/cvt.h"
#include "voronoi/voronoi.h"

using namespace avro;
using namespace avro::delaunay;

UT_TEST_SUITE( cvt_test_suite )

UT_TEST_CASE( test1 )
{
  coord_t number = 3;
  index_t N = 2;

  std::vector<index_t> dims(number,N);
  CKF_Triangulation topology( dims );

  Delaunay sites(topology.points().dim());
  #if 0
  topology.points().copy(sites);
  #else

  // retrieve the number of sites (if random)
  index_t nb_sites = 1e2;

  std::vector<real_t> xmin( sites.dim() ,  1e20 );
  std::vector<real_t> xmax( sites.dim() , -1e20 );
  for (index_t k=0;k<topology.points().nb();k++)
  {
    for (coord_t d=0;d<topology.points().dim();d++)
    {
      if (topology.points()[k][d]<xmin[d]) xmin[d] = topology.points()[k][d];
      if (topology.points()[k][d]>xmax[d]) xmax[d] = topology.points()[k][d];
    }
  }

  std::vector<real_t> x(sites.dim());
  for (index_t k=0;k<nb_sites;k++)
  {
    for (coord_t d=0;d<sites.dim();d++)
      x[d] = random_within( xmin[d] , xmax[d] );
    sites.create(x.data());
  }
  #endif

  delaunay::RestrictedVoronoiDiagram rvd(topology,sites);
  rvd.parallel() = true;

  //rvd.compute(true);
  rvd.optimise(10);

  graphics::Viewer vis;

  //vis.add(topology);
  vis.add(rvd);

  vis.run();
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( cvt_test_suite )
