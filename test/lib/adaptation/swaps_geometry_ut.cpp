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

#include "adaptation/adapt.h"
#include "adaptation/geometry.h"
#include "adaptation/metric.h"
#include "adaptation/parameters.h"

#include "common/error.h"

#include "geometry/egads/context.h"
#include "geometry/tessellation.h"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/egads.h"
#include "library/meshb.h"
#include "library/metric.h"

#include "mesh/mesh.h"

using namespace avro;

UT_TEST_SUITE(swaps_geometry_suite)

UT_TEST_CASE(test1)
{
  // setup the topology
  coord_t number = 2;
  coord_t dim = 3;

  // parameters
  library::MetricField_Uniform analytic(2,0.1);

  // geometry
  EGADS::Context context;
  EGADS::Model model(&context,BASE_TEST_DIR+"/geometry/cube-cylinder.egads");

  TessellationParameters tess_params;
  tess_params.standard();
  tess_params.min_size() = 0.2;
  tess_params.min_angle() = 20;

  ModelTessellation tess(model,tess_params);

  // create a mesh and add the topology
  std::shared_ptr<Mesh> pmesh = std::make_shared<Mesh>(number,dim);
  pmesh->points().set_parameter_dim(number);

  // convert all the points to parameter space
  for (index_t k=0;k<tess.points().nb();k++)
  {
    std::vector<real_t> u = { tess.points().u(k,0) , tess.points().u(k,1) };
    pmesh->points().create( tess.points()[k] );
    pmesh->points().u(k,0) = u[0];
    pmesh->points().u(k,1) = u[1];
    pmesh->points().set_entity(k,tess.points().entity(k));
  }

  // retrieve all the triangles
  Topology<Simplex> topology(pmesh->points(),2);
  tess.retrieve<Simplex>(0).get_elements( topology );
  topology.element().set_parameter(true);

  topology.close();
  topology.neighbours().compute();
  topology.inverse().build();

  EdgeSwap<Simplex> swapper(topology);
  UT_ASSERT( swapper.element().parameter() );

  swapper.delay() = false;

  index_t nb_swaps = 0;
  std::vector<index_t> edges;
  topology.get_edges(edges);
  for (index_t k=0;k<edges.size()/2;k++)
  {
    index_t e0 = edges[2*k];
    index_t e1 = edges[2*k+1];

    std::vector<index_t> elems;
    topology.intersect( {e0,e1},elems);
    UT_ASSERT_EQUALS( elems.size() , 2 );

    // get the entity of this edge
    Entity* entity = swapper.geometry(e0,e1);
    UT_ASSERT( entity!=nullptr );

    // get the opposite vertices
    int pi = topology.neighbours().opposite(elems[0],elems[1]);
    UT_ASSERT( pi >= 0 );
    index_t p = index_t(pi);

    int qi = topology.neighbours().opposite(elems[1],elems[0]);
    UT_ASSERT( qi >= 0 );
    index_t q = index_t(qi);

    bool accept = swapper.valid( p , e0 , e1 );
    if (p==0) UT_ASSERT_EQUALS( accept , false );
    if (!accept)
    {
      continue;
    }

    swapper.set_cavity(elems);

    // apply the swap
    accept = swapper.apply(p,e0,e1);
    if (!accept)
    {
      continue;
    }

    nb_swaps++;

    // swap back, this has to be valid
    accept = swapper.valid( e0 , p , q );
    UT_ASSERT_EQUALS( accept , true );

    /*swapper.set_cavity(elems);
    accept = swapper.apply( e0 , p , q );
    UT_ASSERT_EQUALS( accept , true );*/

    // try the other swap
    accept = swapper.valid( q , e0 , e1 );
    if (p==0) UT_ASSERT_EQUALS( accept , false );
    if (!accept)
    {
      continue;
    }

    // apply the swap
    accept = swapper.apply(p,e0,e1);
    if (!accept)
    {
      continue;
    }

    topology.apply(swapper);

  }

  graphics::Visualizer vis;
  vis.add_topology(topology);

  vis.run();

}
UT_TEST_CASE_END(test1)


UT_TEST_SUITE_END(swaps_geometry_suite)
