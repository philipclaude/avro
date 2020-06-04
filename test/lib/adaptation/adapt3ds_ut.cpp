#include "unit_tester.hpp"

#include "adaptation/adapt.h"
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

UT_TEST_SUITE(adapt3d_surface_suite)

UT_TEST_CASE(test1)
{
  // setup the topology
  coord_t number = 2;
  coord_t dim = 3;

  // parameters
  library::MetricField_Uniform analytic(2,10);

  // geometry
  EGADS::Context context;
  EGADS::Model model(&context,BASE_TEST_DIR+"/geometry/cube-cylinder.egads");

  TessellationParameters tess_params;
  tess_params.standard();
  tess_params.min_size() = 0.5;
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
  std::shared_ptr<Topology<Simplex>> ptopology;
  ptopology = std::make_shared<Topology<Simplex>>(pmesh->points(),2);
  pmesh->add(ptopology);
  tess.retrieve<Simplex>(0).get_elements( *ptopology );
  ptopology->master().set_parameter(true);

  // define the problem and adapt
  AdaptationParameters params;
  params.directory() = "tmp/";
  params.insertion_volume_factor() = -1;
  params.curved() = true;
  params.has_uv() = true;
  params.use_smoothing() = false;
  params.swapout() = false;
  //params.limit_metric() = true; // required a little implementation first

  index_t niter = 2;
  for (index_t iter=0;iter<=niter;iter++)
  {

    params.adapt_iter() = iter;

    // create the metric field
    std::vector<numerics::SymMatrixD<real_t>> fld;
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
    ptopology->master().set_parameter(true);

    if (iter==niter)
    {
      graphics::Visualizer vis;

      //vis.add_topology(topology);
      vis.add_topology(topology_out);

      vis.run();

    }

    params.has_uv() = true;
  }
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(adapt3d_surface_suite)
