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

namespace avro
{
namespace library
{
class MetricField_UGAWG_Linear2 : public MetricField_Analytic
{
public:
  MetricField_UGAWG_Linear2() :
    MetricField_Analytic(2)
  {}

  numerics::SymMatrixD<real_t> operator()( const real_t* x ) const
  {
    numerics::SymMatrixD<real_t> m(dim_);

    real_t hx = 0.1;
    real_t h0 = 1e-3;
    real_t hy = h0 +2.*(0.1 -h0)*fabs( x[1] -0.5 );

    m(0,0) = 1./(hx*hx);
    m(0,1) = 0.;
    m(1,1) = 1./(hy*hy);
    return m;
  }
};
} // library
} // avro

UT_TEST_SUITE( adaptation_parallel_test_suite )

#ifdef AVRO_MPI

UT_TEST_CASE( test1 )
{
  coord_t number = 2;
  coord_t dim = number;

  std::vector<index_t> dims(number,10);
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

  AdaptationParameters params;
  params.standard();

  std::vector<VertexMetric> metrics(topology.points().nb());
  library::MetricField_UGAWG_Linear2 analytic;
  //library::MetricField_Uniform analytic(number,0.2);
  for (index_t k=0;k<topology.points().nb();k++)
    metrics[k] = analytic( topology.points()[k] );

  params.partitioned() = false;
  params.balanced() = true; // assume load-balanced once the first partition is computed
  params.curved() = false;
  params.insertion_volume_factor() = -1;
  params.limit_metric() = true;

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

    // it may not look like much, but this is a huge check on the partitioning and stitching algorithm
    //UT_ASSERT_EQUALS( topology_out.points().nb() , topology.points().nb() );

    graphics::Visualizer vis;
    vis.add_topology(topology_out);
    if (number<4)
      vis.run();
  }

}
UT_TEST_CASE_END( test1 )

#endif

UT_TEST_SUITE_END( adaptation_parallel_test_suite )
