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
class MetricField_UGAWG_sin : public MetricField_Analytic
{
public:
  MetricField_UGAWG_sin() :
    MetricField_Analytic(2)
  {}

  numerics::SymMatrixD<real_t> operator()( const real_t* x ) const
  {
    numerics::SymMatrixD<real_t> m(dim_);

    real_t H_ = 5.;
		real_t f_ = 10.;
		real_t A_ = 2.;
		real_t P_ = 5.;
		real_t y0_ = 5.;
		real_t w_ = 1.*M_PI/P_;

    real_t phi = -f_*( x[1] -A_*cos(w_*x[0]) -y0_ );
  	real_t dzdx = -A_*H_*f_*w_*sin(w_*x[0])*( pow(tanh(phi),2.) -1. );
  	real_t dzdy = -H_*f_*( pow(tanh(phi),2.) -1. );
  	m(0,0) = 1. +dzdx*dzdx;
  	m(0,1) = dzdx*dzdy;
  	m(1,1) = 1. +dzdy*dzdy;
    return m;
  }
};

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
    real_t hy = h0 +2.*(hx -h0)*fabs( x[1] -0.5 );

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
  coord_t number = 3;
  coord_t dim = number;

  std::vector<index_t> dims(number,10);
  CKF_Triangulation topology(dims);

  //for (index_t k=0;k<topology.points().nb();k++)
  //for (coord_t d=0;d<topology.points().dim();d++)
  //  topology.points()[k][d] *= 10;

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
  library::MetricField_UGAWG_Linear analytic;
  //library::MetricField_UGAWG_sin analytic;
  //library::MetricField_Uniform analytic(number,0.2);
  for (index_t k=0;k<topology.points().nb();k++)
    metrics[k] = analytic( topology.points()[k] );

  params.partitioned() = false;
  params.balanced() = true; // assume load-balanced once the first partition is computed
  params.curved() = false;
  params.insertion_volume_factor() = -1;
  params.limit_metric() = true;
  params.max_passes() = 5;
  params.swapout() = false;

  topology.build_structures();

  AdaptationManager<Simplex> manager( topology , metrics , params );

  index_t rank = mpi::rank();

  index_t niter = 5;
  for (index_t iter=0;iter<=niter;iter++)
  {
    if (rank == 0)
      printf("*** global pass %lu ***\n",iter);

    // adapt the mesh, migrate interfaces, etc.
    manager.adapt();

    // re-evaluate the metrics
    metrics.resize( manager.topology().points().nb() );
    for (index_t k=0;k<metrics.size();k++)
      metrics[k] = analytic( manager.topology().points()[k] );

    manager.reassign_metrics(metrics);

    mpi::barrier();
  }

  Points points_out(dim,dim-1);
  Topology<Simplex> topology_out(points_out,number);
  manager.retrieve(topology_out);

  if (rank==0)
  {
    real_t volume = topology_out.volume();
    UT_ASSERT_NEAR( volume , 1.0 , 1e-12 );

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
