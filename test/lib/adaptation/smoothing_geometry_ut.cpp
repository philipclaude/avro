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
  pmesh->points().print(true);

  // retrieve all the triangles
  Topology<Simplex> topology(pmesh->points(),2);
  tess.retrieve<Simplex>(0).get_elements( topology );
  topology.element().set_parameter(true);

  topology.close();
  topology.neighbours().compute();
  topology.inverse().build();

  // create the metric field
  std::vector<symd<real_t>> fld;
  for (index_t k=0;k<topology.points().nb();k++)
    fld.push_back( analytic( topology.points()[k] ) );

  MetricAttachment attachment(topology.points(),fld);
  MetricField<Simplex> metric(topology,attachment);

  Smooth<Simplex> smoother(topology);
  UT_ASSERT( smoother.element().parameter() );

  // loop through every edge, extract the cavity and ensure every edge end point is visible
  for (index_t iter=0;iter<10;iter++)
  for (index_t k=0;k<topology.points().nb();k++)
  {
    if (k < topology.points().nb_ghost()) continue;
    Entity* entity = topology.points().entity(k);
    if (entity->number()!=2) continue;
    bool accept = smoother.apply( k , metric , -1 );
    UNUSED(accept);
  }

  graphics::Visualizer vis;
  vis.add_topology(topology);

  vis.run();

}
UT_TEST_CASE_END(test1)


UT_TEST_SUITE_END(insertions_geometry_suite)
