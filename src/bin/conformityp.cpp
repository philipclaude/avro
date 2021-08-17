#include "adaptation/adapt.h"
#include "adaptation/metric.h"
#include "adaptation/parallel.h"
#include "adaptation/properties.h"

#include "common/mpi.hpp"
#include "common/process.h"
#include "programs.h"

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

namespace avro
{

namespace programs
{

int
conformityp( int nb_input , const char** inputs ) {

  typedef Simplex type;
  avro_assert( nb_input == 6 );

  const char **options = inputs +2;
  int nb_options = nb_input -2;

  // read the base directory and metric name
  std::string base(inputs[0]);
  std::string metric_name(inputs[1]); // such as Polar2, "RotatingBL-4d" , "discrete"
  std::string output(inputs[2]);

  // read the adapt iter
  index_t iter = 5, pass = 2, nb_processors = 1;
  avro_assert( parse<index_t>(lookfor(options,nb_options,"adapt"),iter) );
  avro_assert( parse<index_t>(lookfor(options,nb_options,"pass"),pass) );
  avro_assert( parse<index_t>(lookfor(options,nb_options,"proc"),nb_processors) );

  std::string mesh_name;
  mesh_name = base + "/mesh-adapt" + std::to_string(iter) + "-pass" + std::to_string(pass) + ".avro";

  bool analytic = true;
  real_t scale;

  if (metric_name == "discrete") analytic = false;

  // read in the mesh
  std::shared_ptr<TopologyBase> ptopology = nullptr;
  std::shared_ptr<Mesh> pmesh = library::get_mesh(mesh_name,ptopology);

  // get the topology and add it to the input mesh
  Topology<type>& topology = *static_cast<Topology<type>*>(ptopology.get());
  topology.orient();

  topology.element().set_basis( BasisFunctionCategory_Lagrange );

  index_t n = topology.number();
  index_t nb_rank = n*(n+1)/2;

  // determine the metrics
  std::shared_ptr<MetricAttachment> attachment;
  if (analytic) {

    scale = std::pow( std::sqrt(2.0) , real_t(iter) );

    attachment = library::get_metric( metric_name , topology.points() , analytic );
    avro_assert( analytic );

    // we need to scale the attachment
    for (index_t k = 0; k < attachment->nb(); k++) {

      for (index_t j = 0; j < n; j++)
        for (index_t i = j; i < n; i++)
          (*attachment)[k](i,j) *= scale;
      (*attachment)[k].calculate();
    }
  }
  else {

    scale = 1.0;
    //scale = std::pow( std::sqrt(2.0) , real_t(iter) );


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
          m.data(i) = scale*metrics_p[ j*nb_rank + i ];
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

  properties.dump(output);

  return 0;

}

} // programs

} // avro
