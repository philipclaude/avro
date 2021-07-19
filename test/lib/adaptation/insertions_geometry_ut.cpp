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

UT_TEST_SUITE(insertions_geometry_suite)

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
  tess_params.set("min size" , 0.2 );
  tess_params.set("min angle", 20.0);

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
  pmesh->points().print(true);

  // retrieve all the triangles
  Topology<Simplex> topology(pmesh->points(),2);
  tess.retrieve<Simplex>(0).get_elements( topology );
  topology.element().set_parameter(true);

  topology.close();
  topology.neighbours().compute();
  topology.inverse().build();

  Insert<Simplex> inserter(topology);
  UT_ASSERT( inserter.element().parameter() );

  inserter.delay() = false;

  // loop through every edge, extract the cavity and ensure every edge end point is visible
  std::vector<index_t> edges;
  topology.get_edges(edges);
  for (index_t k=0;k<edges.size()/2;k++)
  {
    index_t p = edges[2*k];
    index_t q = edges[2*k+1];

    // get the entity of this edge
    Entity* entity = inserter.geometry(p,q);
    UT_ASSERT( entity!=nullptr );

    // skip edges along geometry Edges
    //if (entity->number()==1) continue;

    bool accept;

    // compute the midpoint coordinates of an insertion
    std::vector<real_t> U(2);
    std::vector<real_t> X(3);

    // revert the coordinates back to parameter space
    geometry_params( entity , topology.points() , &p , 1 , topology.points()[p] );
    geometry_params( entity , topology.points() , &q , 1 , topology.points()[q] );

    U[0] = 0.5*( topology.points()[p][0] + topology.points()[q][0] );
    U[1] = 0.5*( topology.points()[p][1] + topology.points()[q][1] );
    entity->evaluate(U,X);

    topology.inverse().create(1);

    accept = inserter.apply( p , q , X.data() , U.data() );
    UT_ASSERT( accept) ;
  }

  graphics::Viewer vis;
  vis.add(topology);
  vis.run(AVRO_FULL_UNIT_TEST);

}
UT_TEST_CASE_END(test1)


UT_TEST_SUITE_END(insertions_geometry_suite)
