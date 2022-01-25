#include "unit_tester.hpp"

#include "adaptation/adapt.h"
#include "adaptation/metric.h"
#include "adaptation/parallel.h"
#include "adaptation/properties.h"

#include "common/mpi.hpp"
#include "common/process.h"

#include "geometry/egads/context.h"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/egads.h"
#include "library/factory.h"
#include "library/meshb.h"
#include "library/metric.h"
#include "library/tesseract.h"

#include "mesh/interpolation.h"
#include "mesh/partition.h"

#include "numerics/geometry.h"

#include <fstream>
#include <time.h>

using namespace avro;

UT_TEST_SUITE( adaptation_parallel_test_suite )

#if AVRO_MPI

template<typename F>
void
evaluate_metric( const F& field , const Points& points , std::vector<VertexMetric>& metrics , bool initial_eval ) {
  index_t rank = mpi::rank();
  if (rank == 0)
    printf("evaluating analytic metric\n");
  metrics.resize( points.nb() );
  for (index_t k = 0; k < metrics.size(); k++)
    metrics[k] = field( points[k] );
}

template<>
void
evaluate_metric( const MetricAttachment& attachment , const Points& points , std::vector<VertexMetric>& metrics , bool initial_eval ) {

#if 0
  const ElementSearch<Simplex>& searcher = field.searcher();
  const Field<Simplex,Metric> metric_field = field.field();
  const Topology<Simplex>& topology = metric_field.topology();
  std::vector<real_t> phi( metric_field.element().nb_basis() , 0. );
  std::vector<real_t> xref( topology.element().number() + 1 );

  index_t n = topology.number();
  index_t nb_rank = n*(n+1)/2;

  clock_t TIME0, TIME1;
  real_t search_time = 0.0;

  metrics.resize( points.nb() , VertexMetric(n) );
  for (index_t k = 0; k < metrics.size(); k++) {

    // get the reference coordinates of the closest element
    TIME0 = clock();
    int ielem = searcher.closest( points[k] , xref );
    TIME1 = clock();
    search_time += real_t(TIME1-TIME0)/real_t(CLOCKS_PER_SEC);
    topology.element().reference().basis().evaluate( xref.data() , phi.data() );

    // perform the interpolation and return the element containing the point
    index_t elem = index_t(ielem);

    Metric mk(topology.element().number());
    bool success = metric_field.dof().interpolate( metric_field[elem] , metric_field.nv(elem) , phi , &mk );
    avro_assert( success );

    VertexMetric m(n);
    for (index_t j = 0; j < nb_rank; j++)
      m.data(j) = mk.data(j);
    metrics[k] = m;
  }
  printf("search time = %g\n",search_time);
#else
  if (initial_eval) {
    metrics.resize( points.nb() );
    for (index_t k = 0; k < metrics.size(); k++)
      metrics[k] = attachment[k];
  }
#endif
}

template<typename F>
void
adapt( const F& metric_field , const std::string& mesh_name , const std::string& geometry_name , bool analytic ) {

  typedef Simplex type;

  // get the original input mesh
  std::shared_ptr<TopologyBase> ptopology = nullptr;
  std::shared_ptr<Mesh> pmesh = library::get_mesh(mesh_name,ptopology);

  // get the topology and add it to the input mesh
  Topology<type>& topology = *static_cast<Topology<type>*>(ptopology.get());
  topology.orient();

  // get the input geometry
  bool curved;
  std::shared_ptr<Model> pmodel;
  pmodel = library::get_geometry( geometry_name , curved );
  Model& model = *pmodel;

  // attach the input points to the geometry
  topology.points().attach(model);

  std::vector<VertexMetric> metrics(topology.points().nb());
  evaluate_metric( metric_field , topology.points() , metrics , true );

  AdaptationParameters params;
  params.set( "directory" , std::string("tmp/"));
  params.set( "insertion volume factor" ,  -1.0 );
  params.set( "curved" , curved);
  params.set( "limit metric" , false );
  params.set( "max parallel passes" , index_t(3) );
  params.set( "elems per processor" , index_t(5000) );
  params.set( "has uv", true);
  params.set( "swapout" , false);
  params.set( "force partition count" , index_t(mpi::size()) );
  params.set( "debug level" , index_t(0) );
  params.set( "geometry" , geometry_name );
  params.set( "allow serial" , false );
  params.set( "write conformity" , true );

  topology.build_structures();

  AdaptationManager<Simplex> manager( topology , metrics , params );

  index_t rank = mpi::rank();

  index_t n = topology.number();
  index_t nb_metric_rank = n*(n+1)/2;

  index_t niter = 14;//10;
  for (index_t iter = 0; iter <= niter; iter++) {

    params.set("adapt iter",index_t(iter));

    if (topology.number() == 3 && iter <= 2 && mpi::size() > 16) {
      params.set( "force partition count" , index_t(16) );
    }
    else if (topology.number() == 4 && iter <= 1 && mpi::size() > 4) {
      params.set( "force partition count" , index_t(4) );
    }
    else params.set( "force partition count" , index_t( mpi::size() ) );

    if (rank == 0)
      printf("\n=== iteration %lu ===\n\n",iter);

    // adapt the mesh, migrate interfaces, etc.
    manager.adapt();

    // re-evaluate the metrics
    metrics.resize( manager.topology().points().nb() );
    metrics = manager.metrics();
    real_t scale = sqrt(2.0);
    //if (n < 4) scale = 2.0;
    if (analytic) {
      // we will scale the analytic metric, not the one from the previous adaptation
      scale = std::pow( scale , real_t(iter+1) );
      evaluate_metric( metric_field , manager.topology().points() , metrics , false );
   }

   // save the metric data for this processor
   // the metrics have been migrated so each processor own the metrics for the vertices it needs
   // we need to save the global id for vertex
   std::vector<real_t> metric_data( nb_metric_rank * metrics.size() );
   index_t count = 0;
   for (index_t k = 0; k < metrics.size(); k++)
   for (index_t j = 0; j < nb_metric_rank; j++)
    metric_data[count++] = metrics[k].data(j);

    avro_assert( metrics.size() == manager.topology().points().nb() );
    std::vector<index_t> global( metrics.size() );
    for (index_t k = 0; k < metrics.size(); k++)
      global[k] = manager.topology().points().global(k);

    nlohmann::json jm;
    jm["metrics"] = metric_data;
    jm["global"]  = global;

    std::ofstream file("metric-proc"+std::to_string(rank)+"-iter"+std::to_string(iter)+".metric");
    file << jm;
    file.close();

    // scale the metrics by 2 with respect to the previous iteration
    for (index_t k = 0; k < metrics.size(); k++) {
      for (index_t j = 0; j < metrics[k].m(); j++)
        for (index_t i = j; i < metrics[k].n(); i++)
          metrics[k](i,j) *= scale;
    }
    manager.reassign_metrics(metrics);

    if (rank == 0) {
      printf("==> timing: process = %4.3g sec., adapt = %4.3g sec., synchronize = %4.3g sec.\n",
              manager.time_process(),manager.time_adapt(),manager.time_synchronize());
      printf("            migrate = %4.3g sec., partition = %4.3g sec., exchange = %4.3g sec.\n",
              manager.time_migrate(),manager.time_partition(),manager.time_exchange());

      nlohmann::json jt;
      jt["process"] = manager.time_process();
      jt["adapt"]   = manager.time_adapt();
      jt["synchronize"] = manager.time_synchronize();
      jt["migrate"] = manager.time_migrate();
      jt["partition"] = manager.time_partition();
      jt["exchange"] = manager.time_exchange();

      std::ofstream tfile("timing-iter" + std::to_string(iter) + ".json" );
      tfile << jt;
      tfile.close();
    }
    mpi::barrier();
  }
  fflush(stdout);

}

// which test case to run?
#define CASE_SL 0
#define CASE_CC 0
#define CASE_BL 0
#define CASE_SW 1

UT_TEST_CASE( test1 )
{

  bool analytic = true;
  #if CASE_SL

  std::string geometry_name = "square";
  std::string mesh_name = AVRO_SOURCE_DIR + "/build/release_mpi/sl_9.mesh";
  //std::string mesh_name = "/home/pcaplan/jobs/adapt/ccp2-initial.mesh";
  library::MetricField_UGAWG_Linear2d metric_field;

  #elif CASE_CC
  std::string geometry_name = "box";
  //std::string mesh_name = AVRO_SOURCE_DIR + "/build/release_mpi/ccp2_9.mesh";
  std::string mesh_name = "/home/pcaplan/jobs/adapt/ccp2-initial.mesh";
  library::MetricField_UGAWG_Polar2 metric_field;

  #elif CASE_BL

  std::string geometry_name = "tesseract";
  std::string mesh_name = "/home/pcaplan/jobs/adapt/mbl-initial.avro";
  library::MetricField_Tesseract_RotatingBoundaryLayer metric_field;

  #elif CASE_SW

  analytic = false;
  std::string geometry_name = "tesseract";
  std::string mesh_name = "/home/pcaplan/jobs/adapt/sw-p1-1024.avro";

  // read in the mesh
  std::shared_ptr<TopologyBase> ptopology = nullptr;
  std::shared_ptr<Mesh> pmesh = library::get_mesh(mesh_name,ptopology);
  Topology<Simplex>& topology = *static_cast<Topology<Simplex>*>(ptopology.get());
  index_t n = topology.number();
  index_t nb_rank = n*(n+1)/2;

  // read in the metric
  nlohmann::json jin;
  std::ifstream file_in(mesh_name);
  file_in >> jin;
  std::vector<real_t> data = jin["metric"];

  // create a list of metric tensors
  std::vector<symd<real_t>> metrics( topology.points().nb() );
  avro_assert( data.size() == topology.points().nb()*nb_rank );
  index_t count = 0;
  for (index_t k = 0; k < metrics.size(); k++) {
    symd<real_t> metric(n);
    for (index_t j = 0; j < nb_rank; j++)
      metric.data(j) = data[count++];
    metrics[k] = metric;
  }
  topology.build_structures();
  topology.element().set_basis( BasisFunctionCategory_Lagrange );

  // create the discrete metric field defined on the input mesh
  MetricAttachment metric_field( topology.points() , metrics );
  MetricField<Simplex> field( topology , metric_field );

  Properties properties( topology , field );
  properties.print( "initial metric conformity" );

  #endif

  adapt( metric_field , mesh_name , geometry_name , analytic );

}
UT_TEST_CASE_END( test1 )

#endif

UT_TEST_SUITE_END( adaptation_parallel_test_suite )
