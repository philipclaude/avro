#include "common/error.h"

#include "element/quadrature.h"

#include "mesh/decomposition.h"

#include "numerics/geometry.h"
#include "numerics/integration.h"

#include "voronoi/laguerre.h"
#include "voronoi/voronoi_cell.h"

#include <HLBFGS/HLBFGS.h>

#include <tinymat/types/SurrealS.h>
#include <tinymat/types/SurrealD.h>

#include <nlopt.hpp>

namespace avro
{

namespace delaunay
{

template<typename _T>
class Integrand_Laguerre_Energy : public Integrand<Integrand_Laguerre_Energy<_T>>
{
public:

  typedef _T T;

  Integrand_Laguerre_Energy( const Points& delaunay , const std::vector<index_t>& sites , const std::vector<index_t>& parents ,
    const std::vector<T>& weights, const std::vector<real_t>& mass ) :
    delaunay_(delaunay),
    sites_(sites),
    parents_(parents),
    weights_(weights),
    mass_(mass)
  {}

  bool needs_solution() const { return false; }
  bool needs_gradient() const { return false; }
  bool needs_hessian() const { return false; }

  T operator()( index_t k , const real_t* xref , const real_t* x ) const
  {
    avro_assert_msg( k < parents_.size() , "elem = %lu, nb_parents = %lu", k , parents_.size() );

    index_t site = sites_[parents_[k]];
    avro_assert_msg( site < delaunay_.nb() , "requested site %lu but nb_delaunay = %lu",site,delaunay_.nb());
    const real_t* z = delaunay_[ sites_[parents_[k]] ];
    T f = 0;
    for (coord_t d=0;d<delaunay_.dim();d++)
      f += (z[d] - x[d])*(z[d] - x[d]);
    f -= weights_[ sites_[parents_[k]] ]*weights_[ sites_[parents_[k]] ];
    return f;
  }

  private:
    const Points& delaunay_;
    const std::vector<index_t>& sites_;
    const std::vector<index_t>& parents_;
    const std::vector<T>& weights_;
    const std::vector<real_t>& mass_;
};

LaguerreDiagram::LaguerreDiagram( const Points& sites , const Topology<Simplex>& domain , const std::vector<real_t>& weights ) :
  Topology<Polytope>(points_,domain.number()),
  points_(sites.dim()),
  domain_points_(sites.dim()+1),
  domain_(domain_points_,domain.number()),
  weight_(weights),
  delaunay_(sites.dim()+1),
  exact_(true)
{
  avro_assert( sites.nb() > 0 );
  coord_t dim = sites.dim();

  // copy the domain, considering that it is embedded into dim + 1
  for (index_t k=0;k<domain.points().nb();k++)
  {
    domain_points_.create( domain.points()[k] );
    domain_points_[k][dim] = 0.0;
  }
  domain_.TopologyBase::copy( domain );

  // compute the maximum weight
  if (weight_.size() == 0)
    weight_.resize( sites.nb() , 0. );
  real_t wm = * std::max_element( weight_.begin() , weight_.end() ) + 1e-16;

  for (index_t k=0;k<sites.nb();k++)
  {
    delaunay_.create( sites[k] );
    delaunay_[k][dim] = std::sqrt( wm - weight_[k] );
  }

  mass_.resize( delaunay_.nb() , 0.0 );
}

void
LaguerreDiagram::compute()
{
  // compute the voronoi diagram
  VoronoiDiagram diagram( delaunay_ , domain_ , true );
  diagram.compute(exact_);

  coord_t dim = points_.dim();
  points_.clear();
  points_.set_dim(dim);
  clear();

  // copy the points, removing the extra coordinate (in dim+1)
  for (index_t k=0;k<diagram.points().nb();k++)
  {
    points_.create( diagram.points()[k] );
    std::vector<int> bisectors = diagram.points().incidence().get(k);
    points_.incidence().add( bisectors.data() , bisectors.size() );
  }
  TopologyBase::copy( diagram );
  sites_ = diagram.sites();

  avro_assert( sites_.size() == nb() );
}

real_t
LaguerreDiagram::eval_objective( std::vector<real_t>& dE_dZ , std::vector<real_t>& dE_dW ) const
{
  // initialize size of centroids
  coord_t dim = points_.dim();
  Points centroids( dim );
  std::vector<real_t> x0( dim , 0. );
  for (index_t k=0;k<delaunay_.nb();k++)
    centroids.create( x0.data() );

  std::vector<real_t> V( delaunay_.nb() , 0. );

  SimplicialDecomposition<Polytope> decomposition(*this);
  decomposition.extract();

  std::vector<index_t> simplices;
  std::vector<index_t> parents;
  decomposition.get_simplices( number() , simplices , parents );
  index_t nb_simplices = simplices.size()/(number_+1);

  Topology<Simplex> simplex_topology( decomposition.points() , number_ );
  const Simplex& element = decomposition.element();
  for (index_t k=0;k<nb_simplices;k++)
  {
    simplex_topology.add( &simplices[k*(number_+1)] , number_+1 );

    real_t vk = element.volume( decomposition.points() , &simplices[k*(number_+1)] , number_+1 );

    index_t s = sites_[ parents[k] ];
    avro_assert( s < V.size() );
    V[s] += vk;

    // compute the centroid of this piece
    std::vector<real_t> c( dim );
    numerics::centroid( &simplices[k*(number_+1)] , number_+1 , decomposition.points() , c );

    // add the contribution to x*V and V
    for (coord_t d=0;d<delaunay_.dim();d++)
      centroids[s][d] += c[d]*vk;
  }

  // go back and divide by the volume
  for (index_t k=0;k<centroids.nb();k++)
  for (coord_t d=0;d<centroids.dim();d++)
  {
    // set the centroid to the original seed if the cell volume is zero
    if (V[k]<1e-12)
      centroids[k][d] = delaunay_[k][d];
    else
      centroids[k][d] /= V[k];
  }

  // calculate the gradients of the energy w.r.t. delaunay points and weights
  dE_dZ.resize( delaunay_.nb()*dim , 0. );
  dE_dW.resize( delaunay_.nb() , 0. );

  real_t gnorm_x = 0.;
  real_t gnorm_w = 0.;
  for (index_t k=0;k<delaunay_.nb();k++)
  {
    for (index_t d=0;d<dim;d++)
    {
      dE_dZ[ k*dim + d ] = 2.0*V[k]*( delaunay_[k][d] - centroids[k][d] );
      gnorm_x += std::pow(dE_dZ[k*dim+d],2);
    }

    dE_dW[k] = -V[k] + mass_[k];
    gnorm_w  += std::pow(dE_dW[k],2);
  }

  gnorm_x = std::sqrt( gnorm_x );
  gnorm_w = std::sqrt( gnorm_w );

  // compute the CVT energy
  ConicalProductQuadrature quadrature(number_,3);
  simplex_topology.element().set_basis( BasisFunctionCategory_Lagrange );
  quadrature.define();
  simplex_topology.element().load_quadrature(quadrature);

  typedef Integrand_Laguerre_Energy<real_t> Integrand_t;
  Integrand_t integrand(delaunay_,sites_,parents,weight_,mass_);

  Functional<Integrand_t> f(integrand);
  f.integrate( simplex_topology );

  printf("f = %g, |g_x| = %g, |g_w| = %g\n",f.value(),gnorm_x,gnorm_w);

  real_t energy = f.value();
  for (index_t k=0;k<delaunay_.nb();k++)
    energy += weight_[k]*mass_[k];

  return energy;
}

static void
hlbfgs_eval_func( int n , double* x , double* prev_x , double *f , double* g )
{

}

static void
hlbfgs_new_iteration( int iter , int call_iter , double* x , double* f , double* g , double* gnorm )
{
  UNUSED(iter);
  UNUSED(call_iter);
  UNUSED(x);
  UNUSED(f);
  UNUSED(g);
  UNUSED(gnorm);
}

struct nlopt_data_cvt
{
	LaguerreDiagram& laguerre;
	index_t eval_count;
	real_t objective;
};

real_t
LaguerreDiagram::eval_objective() const
{
  std::vector<real_t> dE_dZ,dE_dW;
  return eval_objective(dE_dZ,dE_dW);
}

void
LaguerreDiagram::set_delaunay( const real_t* x , coord_t dim )
{
  for (index_t k=0;k<delaunay_.nb();k++)
  for (coord_t d=0;d<dim;d++)
    delaunay_[k][d] = x[k*dim+d];
}

double
nlopt_cvt_objective( unsigned n , const double* x , double* grad, void* data0 )
{
  nlopt_data_cvt* data = (nlopt_data_cvt*)(data0);

  LaguerreDiagram& laguerre = data->laguerre;

  // set the coordinates of the delaunay points
  laguerre.set_delaunay(x,laguerre.points().dim());

  // recompute the voronoi diagram
  laguerre.compute();

  // compute the energy and gradients
  std::vector<real_t> dE_dZ,dE_dW;
  real_t energy = laguerre.eval_objective( dE_dZ , dE_dW );
  if (grad != nullptr)
  {
    for (index_t k=0;k<dE_dZ.size();k++)
      grad[k] = dE_dZ[k];
  }

  return energy;
}

void
LaguerreDiagram::optimize_cvt()
{
  coord_t dim = points_.dim();
  index_t n = delaunay_.nb()*dim;

  std::vector<real_t> x(n);
  for (index_t k=0;k<delaunay_.nb();k++)
  for (index_t d=0;d<dim;d++)
    x[k*dim+d] = delaunay_(k,d);

  // setup the optimizer
  nlopt::opt opt( nlopt::LD_LBFGS , n );
  nlopt_data_cvt data = {*this,0,1};

	// set the objective function
	opt.set_min_objective( &nlopt_cvt_objective , static_cast<void*>(&data) );

  // set some optimization parameters
  opt.set_stopval(1e-7);
  opt.set_xtol_rel(1e-7);
  opt.set_ftol_rel(1e-7);
  opt.set_maxeval(20);

  opt.optimize(x);

}

void
LaguerreDiagram::optimize_otm()
{

}

} // delaunay

} // avro
