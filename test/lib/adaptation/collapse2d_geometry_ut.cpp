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

#include "adaptation/adapt.h"
#include "adaptation/geometry.h"
#include "adaptation/metric.h"

#include "common/error.h"

#include "geometry/egads/context.h"
#include "geometry/tessellation.h"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/egads.h"
#include "library/factory.h"
#include "library/meshb.h"
#include "library/metric.h"

#include "mesh/mesh.h"

using namespace avro;

UT_TEST_SUITE(collapses_geometry_suite)

UT_TEST_CASE(test1)
{
  typedef Simplex type;

  // setup the topology
  coord_t number = 2;
  coord_t dim = 2;

  const std::string meshname = AVRO_SOURCE_DIR + "/test/library/meshes/cyl2d.meshb";
  const std::string geometryname = AVRO_SOURCE_DIR + "/test/library/geometry/cyl2d.egads";

  // get the original input mesh
  std::shared_ptr<TopologyBase> ptopology = nullptr;
  std::shared_ptr<Mesh> pmesh = library::get_mesh(meshname,ptopology);
  UT_ASSERT_EQUALS( pmesh->number() , number );

  // get the topology and add it to the input mesh
  Topology<type>& topology = *static_cast<Topology<type>*>(ptopology.get());
  topology.orient();

  // get the input geometry
  bool curved = false;
  std::shared_ptr<Model> pmodel;
  pmodel = library::get_geometry( geometryname , curved );
  Model& model = *pmodel;
  UT_ASSERT( curved );

  // check the points are on the geometry...
  // option to project them
  pmesh->points().attach(model);

  // setup the topology
  topology.close();
  topology.neighbours().compute();
  topology.inverse().build();

  for (index_t k = 0; k < topology.points().nb(); k++) {
    Entity* e = topology.points().entity(k);
    if (e == nullptr) continue;
    //if (e->number() == 2) topology.points().set_entity(k,nullptr);
  }

  topology.points().print(true);

  Collapse<Simplex> collapser(topology);

  collapser.delay() = false;

  index_t pass = 0;
  while (true) {

    pass++;
    std::vector<bool> removed( topology.points().nb() , false );

    index_t nb_collapse = 0;

    // loop through every edge, extract the cavity and ensure every edge end point is visible
    std::vector<index_t> edges;
    topology.get_edges(edges);
    for (index_t k = 0; k < edges.size()/2; k++) {

      index_t p0 = edges[2*k];
      index_t p1 = edges[2*k+1];

      if (removed[p0] || removed[p1]) continue;

      // skip collapses with ghost points
      if (p0 < topology.points().nb_ghost() || p1 < topology.points().nb_ghost())
        continue;

      // skip collapses with ghost points
      if (p0 < topology.points().nb_ghost() || p1 < topology.points().nb_ghost())
        continue;

      // get the entity of this edge
      //Entity* entity = collapser.geometry(p0,p1);

      // check if the collapse is valid
      bool accept = false;
      if (collapser.valid(p0,p1)) {
        accept = collapser.apply( p0 , p1 );
        if (accept) removed[p0] = true;
      }
      else if (collapser.valid(p1,p0)) {
        accept = collapser.apply( p1 , p0 );
        if (accept) removed[p1] = true;
      }
      if (accept) nb_collapse++;
    }
    if (nb_collapse == 0) break;
  }

  #if 0
  // now let's try inserting!
  Insert<Simplex> inserter(topology);
  for (index_t pass = 0; pass < 5; pass++) {

    index_t nb_insert = 0;

    // loop through every edge, extract the cavity and ensure every edge end point is visible
    std::vector<index_t> edges;
    topology.get_edges(edges);
    for (index_t k = 0; k < edges.size()/2; k++) {

      index_t p0 = edges[2*k];
      index_t p1 = edges[2*k+1];

      // skip collapses with ghost points
      if (p0 < topology.points().nb_ghost() || p1 < topology.points().nb_ghost())
        continue;

      // skip collapses with ghost points
      if (p0 < topology.points().nb_ghost() || p1 < topology.points().nb_ghost())
        continue;

      // get the entity of this edge
      Entity* entity = collapser.geometry(p0,p1);
      if (entity != nullptr && entity->number() != 1) continue;

      real_t u0 = topology.points().u(p0)[0];
      real_t u1 = topology.points().u(p1)[0];
      std::vector<real_t> u = {0.5*(u0 + u1)};

      // insert on the edge
      std::vector<real_t> X(dim);
      entity->evaluate(u,X);

      topology.inverse().create(1);

      bool accept = inserter.apply( p0 , p1 , X.data() , u.data() );
      UT_ASSERT( accept) ;
    }
  }
  #endif

  graphics::Viewer vis;
  vis.add(topology);
  vis.run(AVRO_FULL_UNIT_TEST);

}
UT_TEST_CASE_END(test1)


UT_TEST_SUITE_END(collapses_geometry_suite)
