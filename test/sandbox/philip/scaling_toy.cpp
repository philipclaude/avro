#include "unit_tester.hpp"

#include "adaptation/adapt.h"
#include "adaptation/parallel.h"

#include "common/mpi.hpp"
#include "common/process.h"

#include "geometry/egads/context.h"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/egads.h"
#include "library/meshb.h"
#include "library/metric.h"
#include "library/tesseract.h"

#include "mesh/partition.h"

#include "numerics/geometry.h"

using namespace avro;

UT_TEST_SUITE( adaptation_parallel_test_suite )

#if AVRO_MPI

UT_TEST_CASE( test1 )
{
  coord_t number = 0, dim = 0;

  EGADS::Context context;

  std::string geometry_name;
  std::string mesh_name;
  #if 1
  dim = number = 3;
  geometry_name = "box";
  mesh_name = AVRO_SOURCE_DIR + "/build/release_mpi/ccp2-box_9.mesh";
  library::MetricField_UGAWG_Polar2 analytic;
  #elif 0
  dim = number = 2;
  geometry_name = "square";
  mesh_name = AVRO_SOURCE_DIR + "/build/release_mpi/sl_19.mesh";
  library::MetricField_UGAWG_Linear2d analytic;
  #else
  #error "unknown test case"
  #endif

  std::vector<real_t> lens(number,1.0);
  EGADS::Cube geometry(&context,lens);
  library::meshb mesh(mesh_name);
  std::shared_ptr<Topology<Simplex>> ptopology = mesh.retrieve_ptr<Simplex>(0);
  Topology<Simplex>& topology = *ptopology.get();

  avro_assert_msg( number > 0 && dim > 0 , "bad dimension" );

  topology.points().attach(geometry);

  std::vector<VertexMetric> metrics(topology.points().nb());
  for (index_t k = 0; k < topology.points().nb(); k++)
    metrics[k] = analytic( topology.points()[k] );

  AdaptationParameters params;
  params.set( "directory" , std::string("tmp/"));
  params.set( "insertion volume factor" ,  -1.0 );
  params.set( "curved" , false);
  params.set( "limit metric" , false );
  params.set( "max parallel passes" , index_t(3) );
  params.set( "elems per processor" , index_t(5000) );
  params.set( "has uv", true);
  params.set( "swapout" , false);
  params.set( "force partition count" , index_t(mpi::size()) );
  params.set( "debug level" , index_t(0) );
  params.set( "geometry" , geometry_name );
  params.set( "allow serial" , false );

  topology.build_structures();

  AdaptationManager<Simplex> manager( topology , metrics , params );

  index_t rank = mpi::rank();

  index_t niter = 10;
  for (index_t iter = 0; iter <= niter; iter++) {

    params.set("adapt iter",index_t(iter));

    if (rank == 0)
      printf("\n=== iteration %lu ===\n\n",iter);

    // adapt the mesh, migrate interfaces, etc.
    manager.adapt();

    // re-evaluate the metrics
    metrics.resize( manager.topology().points().nb() );
    for (index_t k = 0; k < metrics.size(); k++)
      metrics[k] = analytic( manager.topology().points()[k] );

    // scale the metrics by 2 with respect to the previous iteration
    for (index_t k = 0; k < metrics.size(); k++) {
      for (index_t j = 0; j < metrics[k].m(); j++)
        for (index_t i = j; i < metrics[k].n(); i++)
          metrics[k](i,j) *= std::pow(2.0,real_t(iter));
    }

    manager.reassign_metrics(metrics);

    if (rank == 0) {
    printf("==> timing: process = %4.3g sec., adapt = %4.3g sec., synchronize = %4.3g sec.\n",
            manager.time_process(),manager.time_adapt(),manager.time_synchronize());
    printf("            migrate = %4.3g sec., partition = %4.3g sec., exchange = %4.3g sec.\n",
            manager.time_migrate(),manager.time_partition(),manager.time_exchange());
    }
    mpi::barrier();
  }
  fflush(stdout);

}
UT_TEST_CASE_END( test1 )

#endif

UT_TEST_SUITE_END( adaptation_parallel_test_suite )
