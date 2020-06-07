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

#include "common/tools.h"

#include "geometry/entity.h"
#include "geometry/egads/context.h"
#include "geometry/egads/object.h"

#include "graphics/application.h"
#include "graphics/user_interface.h"
#include "graphics/window.h"

#include "library/ckf.h"
#include "library/egads.h"
#include "library/plots.h"

#include "mesh/decomposition.h"

#include "voronoi/delaunay.h"
#include "voronoi/geometry.h"
#include "voronoi/voronoi.h"

using namespace avro;
using namespace avro::graphics;

UT_TEST_SUITE( geometry_rvd_test_suite )

UT_TEST_CASE( test_2d )
{
  coord_t number = 2;
  coord_t dim = number;
  index_t N = 5;
  std::vector<index_t> dims(number,N);
  CKF_Triangulation topology( dims );

  EGADS::Context context;
#if 1
  EGADS::Cube geometry(&context,{1,1});
#else
  real_t xc[3] = {0.5,0.5,0.0};
  EGADS::Square geometry(&context,xc,1.,1.);
#endif
  topology.points().attach(geometry);

  topology.points().print(true);

  Delaunay delaunay(topology.points().dim());
  for (index_t k=0;k<topology.points().nb();k++)
  {
    if (topology.points().entity(k)==nullptr) continue;
    delaunay.create( topology.points()[k] );
    delaunay.set_entity( delaunay.nb()-1 , topology.points().entity(k) );
  }

  index_t ns = 100;
  for (index_t k=0;k<ns;k++)
  {
    std::vector<real_t> x(dim);
    for (coord_t d=0;d<dim;d++)
      x[d] = random_within(0.01,0.99);
    delaunay.create(x.data());
  }

  //topology.points().copy(delaunay);

  delaunay.print(true);

  GeometryConformingRVD geometry_rvd( topology , delaunay );
  geometry_rvd.compute();
  geometry_rvd.extract_triangulations();

  library::Plot<Simplex> point_plot(delaunay);

  Visualizer vis;
  std::shared_ptr<Widget> toolbar = std::make_shared<Toolbar>(vis.main_window(),vis);
  vis.main_window().interface().add_widget( toolbar );

  vis.add_topology(geometry_rvd);
  vis.add_topology(geometry_rvd.triangulation());
  vis.add_topology(point_plot);
  vis.run();
}
UT_TEST_CASE_END( test_2d )

UT_TEST_SUITE_END( geometry_rvd_test_suite )
