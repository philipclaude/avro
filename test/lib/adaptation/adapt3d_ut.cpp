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
#include "adaptation/metric.h"

#include "common/error.h"

#include "geometry/egads/context.h"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/egads.h"
#include "library/meshb.h"
#include "library/metric.h"

#include "mesh/mesh.h"

using namespace avro;

UT_TEST_SUITE( adaptation_adapt3d_suite)

UT_TEST_CASE(adapt_test)
{
  // setup the topology
  coord_t number = 3;

  EGADS::Context context;
  std::vector<real_t> lens(number,1.0);
  std::vector<real_t> x0(number,0.0);

  // parameters
  library::MetricField_UGAWG_Linear analytic;

  // geometry
  std::vector<real_t> lengths(number,1);
  //EGADS::Cube box(&context,lengths);

  EGADS::Model model(&context,BASE_TEST_DIR+"/geometry/cube-cylinder.egads");
  Body& box = model.body(0);

  // structured grid
  std::vector<index_t> dims(number,6);
  std::shared_ptr<Topology<Simplex>> ptopology;
  //ptopology = std::make_shared<CKF_Triangulation>(dims);

  library::meshb mesh(BASE_TEST_DIR+"/meshes/cube-cylinder.mesh");
  ptopology = mesh.retrieve_ptr<Simplex>(0);

  // tag the points onto the body
  ptopology->points().attach( box );

  // create a mesh and add the topology
  std::shared_ptr<Mesh> pmesh = std::make_shared<Mesh>(number,number);
  pmesh->add(ptopology);
  ptopology->points().copy(pmesh->points());

  // define the problem and adapt
  AdaptationParameters params;
  params.set( "directory" , std::string("tmp/"));
  params.set( "insertion volume factor" ,  -1.0 );
  params.set( "curved" , true);
  params.set( "limit metric" , true );

  index_t niter = 1;
  for (index_t iter = 0; iter < niter; iter++)
  {

    params.set( "adapt iter" , index_t(iter) );

    // create the metric field
    std::vector<symd<real_t>> fld;
    for (index_t k=0;k<pmesh->points().nb();k++)
      fld.push_back( analytic( pmesh->points()[k] ) );

    // create the mesh we will write to
    std::shared_ptr<Mesh> pmesh_out = std::make_shared<Mesh>(number,number);

    // define the problem and adapt
    AdaptationProblem problem = {*pmesh,fld,params,*pmesh_out};
    adapt<Simplex>( problem );

    // create the mesh for the next adaptation
    pmesh = std::make_shared<Mesh>(number,number);
    pmesh_out->points().copy( pmesh->points() );
    UT_ASSERT_EQUALS( pmesh->points().nb_ghost() , 0 );

    std::shared_ptr<Topology<Simplex>> ptopology = std::make_shared<Topology<Simplex>>(pmesh->points(),number);
    const Topology<Simplex>& topology_out = pmesh_out->retrieve<Simplex>(0);
    for (index_t k=0;k<topology_out.nb();k++)
      ptopology->add(topology_out(k),topology_out.nv(k));
    pmesh->add(ptopology);

    if (iter==niter)
    {
      graphics::Visualizer vis;

      //vis.add_topology(topology);
      vis.add_topology(topology_out);

      //vis.run();

    }

    params.set("has uv", true);
  }
}
UT_TEST_CASE_END(adapt_test)

UT_TEST_SUITE_END(adaptation_adapt3d_suite)
