#include "unit_tester.hpp"

#include "adaptation/adapt.h"
#include "adaptation/parameters.h"

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

#include "voronoi/delaunay.h"
#include "measures.h"

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

    real_t hu = 0.1;
    real_t h0 = hu/100;
    real_t hy = h0 +2.*(hu -h0)*fabs( x[1] -0.5 );
    real_t hx = hu;//h0 +2.*(hu -h0)*fabs( x[0] -0.5 );

    hx = hy;

    m(0,0) = 1./(hx*hx);
    m(0,1) = 0.;
    m(1,1) = 1./(hy*hy);
    return m;
  }
};

class MetricField_exp : public MetricField_Analytic
{
public:
  MetricField_exp() :
    MetricField_Analytic(2)
  {}

  numerics::SymMatrixD<real_t> operator()( const real_t* x ) const
  {
    numerics::SymMatrixD<real_t> m(dim_);

    real_t c[2] = {0.5,0.5};
    real_t r = std::sqrt( (x[0]-c[0])*(x[0]-c[0]) + (x[1]-c[1])*(x[1]-c[1]) );

    real_t a = 0.05;
    real_t h = a * std::exp( -1 * r*r );

    m(0,0) = 1./(h*h);
    m(0,1) = 0.;
    m(1,1) = 1./(h*h);
    return m;
  }
};

class MetricField_Gaussian : public MetricField_Analytic
{
public:
  MetricField_Gaussian( const tinymat::DLA::VectorD<real_t>& mu , const numerics::SymMatrixD<real_t>& sigma ) :
    MetricField_Analytic(mu.n()),
    density_(mu,sigma)
  {}

  numerics::SymMatrixD<real_t> operator()( const real_t* x ) const
  {
    numerics::SymMatrixD<real_t> m(2);
    real_t rho = 5*density_.evaluate(0,nullptr,x);
    m(0,0) = rho*rho;
    m(0,1) = 0.;
    m(1,1) = rho*rho;
    return m;
  }

private:
  const DensityMeasure_Gaussian density_;
};

} // library
} // avro

UT_TEST_SUITE( adaptation_parallel_test_suite )

UT_TEST_CASE( test1 )
{
  coord_t number = 2;

  EGADS::Context context;
  std::vector<real_t> lengths(number,1.0);
  #if 1
  EGADS::Cube geometry(&context,lengths);

  // gaussian
  tinymat::DLA::VectorD<real_t> mu(number);
  numerics::SymMatrixD<real_t> sigma(number,number);
  sigma = 0;
  for (coord_t d = 0; d < number; d++)
  {
    mu(d) = 0.5;
    sigma(d,d) = 0.05;
  }
  library::MetricField_Gaussian analytic(mu,sigma);
  #else
  std::vector<real_t> c(4,0.5);
  library::Tesseract geometry(c,lengths);
  library::MetricField_Tesseract_Linear analytic;
  #endif

  std::vector<index_t> dims(number,10);
  std::shared_ptr<Topology<Simplex>> ptopology = std::make_shared<CKF_Triangulation>(dims);
  ptopology->points().attach( geometry );

  AdaptationParameters params;
  params.standard();

  params.directory() = "tmp/";
  params.curved() = false;
  params.insertion_volume_factor() = -1;
  params.limit_metric() = true;
  params.swapout() = true;
  params.has_uv() = false;
  params.use_smoothing() = true;

  // create a mesh and add the topology
  std::shared_ptr<Mesh> pmesh = std::make_shared<Mesh>(number,number);
  pmesh->add(ptopology);
  ptopology->points().copy(pmesh->points());

  index_t niter = 5;
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
UT_TEST_CASE_END( test1 )


UT_TEST_SUITE_END( adaptation_parallel_test_suite )
