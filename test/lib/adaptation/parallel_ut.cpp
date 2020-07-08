#include "unit_tester.hpp"

#include "adaptation/parallel.h"
#include "adaptation/parameters.h"

#include "common/mpi.hpp"
#include "common/process.h"

#include "geometry/egads/context.h"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/egads.h"
#include "library/metric.h"
#include "library/tesseract.h"

#include "mesh/partition.h"

#include "numerics/geometry.h"

using namespace avro;

UT_TEST_SUITE( adaptation_parallel_test_suite )

#ifdef AVRO_MPI

UT_TEST_CASE( test1 )
{
  coord_t number = 3;
  coord_t dim = number;

  std::vector<index_t> dims(number,6);
  CKF_Triangulation topology(dims);

  #if 1
  EGADS::Context context;
  std::vector<real_t> lens(number,1.);
  EGADS::Cube geometry(&context,lens);
  #else
  std::vector<real_t> c(4,0.5);
  std::vector<real_t> lengths(4,1.0);
  library::Tesseract geometry(c,lengths);
  #endif
  topology.points().attach(geometry);

  topology.neighbours().fromscratch() = true;
  topology.neighbours().compute();

  AdaptationParameters params;
  params.standard();

  std::vector<VertexMetric> metrics(topology.points().nb());
  library::MetricField_UGAWG_Linear analytic;
  for (index_t k=0;k<topology.points().nb();k++)
    metrics[k] = analytic( topology.points()[k] );

  params.partitioned() = false;
  params.balanced() = true; // assume load-balanced once the first partition is computed
  params.curved() = false;
  params.insertion_volume_factor() = -1;
  //params.limit_metric() = false;

  topology.build_structures();

  AdaptationManager<Simplex> manager( topology , metrics , params );
  manager.adapt();

  Points points_out(dim,dim-1);
  Topology<Simplex> topology_out(points_out,number);
  manager.retrieve(topology_out);

  index_t rank = mpi::rank();
  if (rank==0)
  {
    real_t volume = topology_out.volume();
    //UT_ASSERT_NEAR( volume , 1.0 , 1e-12 );

    for (index_t k=0;k<topology_out.points().nb();k++)
    for (index_t j=k+1;j<topology_out.points().nb();j++)
    {
      real_t d = numerics::distance( topology_out.points()[k] , topology_out.points()[j] , dim );
      if (d < 1e-12)
      {
        //topology_out.points().print(k,true);
        //topology_out.points().print(j,true);
      }
    }

    // it may not look like much, but this is a huge check on the partitioning and stitching algorithm
    //UT_ASSERT_EQUALS( topology_out.points().nb() , topology.points().nb() );

    printf("original |points| = %lu, new |points| = %lu\n",topology.points().nb(),topology_out.points().nb());

    graphics::Visualizer vis;
    vis.add_topology(topology_out);
    if (number<4)
      vis.run();
  }

}
UT_TEST_CASE_END( test1 )

#endif

UT_TEST_SUITE_END( adaptation_parallel_test_suite )
