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

UT_TEST_CASE( test1 )
{

  bool analytic = true;
  index_t iter = 2;
  index_t nb_processors = 24;
  real_t scale;

  std::string mesh_name;
  std::string metric_name;

  mesh_name = "/home/pcaplan/jobs/adapt/ccp2/mesh-adapt" + std::to_string(iter)+ "-pass1.avro";
  metric_name = "Polar2";
 
  mesh_name = "/home/pcaplan/jobs/adapt/results/bl-p24/mesh-adapt" + std::to_string(iter) + "-pass1.avro"; 
  //metric_name = "RotatingBL-4d"; 

  std::string base = "/home/pcaplan/jobs/adapt/results/bl-p24/";
  //mesh_name = "/home/pcaplan/jobs/adapt/sw/mesh-adapt" + std::to_string(iter) + "-pass1.avro";
  metric_name = "discrete";

  typedef Simplex type;

  // read in the mesh
  std::shared_ptr<TopologyBase> ptopology = nullptr;
  std::shared_ptr<Mesh> pmesh = library::get_mesh(mesh_name,ptopology);

  // get the topology and add it to the input mesh
  Topology<type>& topology = *static_cast<Topology<type>*>(ptopology.get());
  topology.orient();

  topology.build_structures();
  topology.element().set_basis( BasisFunctionCategory_Lagrange );
  
  index_t n = topology.number();
  index_t nb_rank = n*(n+1)/2;


  // determine the metrics
  std::shared_ptr<MetricAttachment> attachment;
  try {
    attachment = library::get_metric( metric_name , topology.points() , analytic );
    avro_assert( analytic );
  
    scale = std::pow( std::sqrt(2.0) , real_t(iter) );

    // we need to scale the attachment
    for (index_t k = 0; k < attachment->nb(); k++) {

      for (index_t j = 0; j < n; j++)
        for (index_t i = j; i < n; i++)
          (*attachment)[k](i,j) *= scale;
      (*attachment)[k].calculate();
    }

  }
  catch (...) {

   // we need to read the metrics for each processor
   std::vector< symd<real_t> > metrics( topology.points().nb() );
   std::unordered_set<index_t> added; 
   for (index_t k = 0; k < nb_processors; k++) {

	std::string metric_name_k = base + "/metric-proc" + std::to_string(k) + "-iter" + std::to_string(iter) + ".metric";

        printf("reading metric %s from processor %lu\n",metric_name_k.c_str(),k);

	nlohmann::json jm;
        std::ifstream file_in(metric_name_k);
        file_in >> jm;

        std::vector<real_t> metrics_p = jm["metrics"];
        std::vector<real_t> global_p  = jm["global"];
  
        avro_assert( global_p.size() * nb_rank == metrics_p.size() );

        // append the metrics into the global list
        for (index_t j = 0; j < global_p.size(); j++) {

          index_t idx = global_p[j];
          if (added.find(idx) != added.end()) continue;
          added.insert(idx);
 
          symd<real_t> m(n);
          for (index_t i = 0; i < nb_rank; i++)
            m.data(i) = metrics_p[ j*nb_rank + i ];
          metrics[idx] = m;
        }
    }
    printf("done reading metrics\n");
    avro_assert( added.size() == topology.points().nb() );
    attachment = std::make_shared<MetricAttachment>( topology.points() , metrics );
  }


  // create the discrete metric field defined on the input mesh
  MetricField<Simplex> field( topology , *attachment );

  Properties properties( topology , field );
  properties.print( "metric conformity" );

}
UT_TEST_CASE_END( test1 )

#endif

UT_TEST_SUITE_END( adaptation_parallel_test_suite )
