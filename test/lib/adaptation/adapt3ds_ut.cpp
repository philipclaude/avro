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
#include "geometry/tessellation.h"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/csm.h"
#include "library/egads.h"
#include "library/meshb.h"
#include "library/metric.h"

#include "mesh/mesh.h"

using namespace avro;

UT_TEST_SUITE(adapt3d_surface_suite)

UT_TEST_CASE(test1)
{
  // setup the topology
  coord_t number = 2;
  coord_t dim = 3;

  // parameters
  std::shared_ptr<library::MetricField_UniformGeometry<Simplex>> metric;
  metric = std::make_shared<library::MetricField_UniformGeometry<Simplex>>(2,0.1);

  // geometry
  #if 1
  EGADS::Context context;
  EGADS::Model model(&context,BASE_TEST_DIR+"/geometry/cube-cylinder.egads");
  #elif 0 // doesn't work because I think Face3D is negatively oriented
  EGADS::Context context;
  std::shared_ptr<Body> face = std::make_shared<EGADS::Face3D>( &context , 0 );
  EGADS::Model model(2);
  model.add_body(face);
  #elif 1
  EGADS::Context context;
  real_t x0[3] = {0,0,0};
  real_t dir[3] = {0,0,1};
  std::shared_ptr<Body> face = std::make_shared<EGADS::Torus>( &context , x0,dir,1.,0.1);
  EGADS::Model model(2);
  model.add_body(face);
  #else
  OpenCSM_Model model("/Users/pcaplan/Codes/EngSketchPad/data/bottle.csm");
  #endif

  TessellationParameters tess_params;
  tess_params.set("min size" , 1.0 );
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
  //pmesh->points().print(true);

  // retrieve all the triangles
  std::shared_ptr<Topology<Simplex>> ptopology;
  ptopology = std::make_shared<Topology<Simplex>>(pmesh->points(),2);
  pmesh->add(ptopology);
  tess.retrieve<Simplex>(0).get_elements( *ptopology );
  ptopology->element().set_parameter(true);
  ptopology->element().set_basis( BasisFunctionCategory_Lagrange );

  // define the problem and adapt
  AdaptationParameters params;
  params.set( "directory" , std::string("tmp/"));
  params.set( "insertion volume factor" ,  -1.0 );
  params.set( "curved" , true);
  params.set( "limit metric" , false ); // requires some implementation
  params.set( "has uv",  true );
  params.set( "use smoothing" , true);
  params.set( "swapout" , false);

  index_t niter = 1;
  for (index_t iter = 0; iter < niter; iter++)
  {

    params.set( "adapt iter" , index_t(iter) );

    // create the metric field
    std::vector<symd<real_t>> fld;
//    #pragma omp parallel for
    for (index_t k=0;k<pmesh->points().nb();k++)
      fld.push_back( (*metric.get())( pmesh->points() , k ) );

    // create the mesh we will write to
    std::shared_ptr<Mesh> pmesh_out = std::make_shared<Mesh>(number,number);

    // define the problem and adapt
    AdaptationProblem problem = {*pmesh,fld,params,*pmesh_out , metric.get() };
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
    ptopology->element().set_parameter(true);

    if (iter==niter)
    {
      graphics::Visualizer vis;
      vis.add_topology(topology_out);
      //vis.run();
    }

    params.set("has uv", true);
  }
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(adapt3d_surface_suite)
