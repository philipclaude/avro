#include "common/error.h"
#include "common/process.h"

#include "element/quadrature.h"

#include "mesh/decomposition.h"

#include "numerics/geometry.h"
#include "numerics/integration.h"
#include "numerics/nlopt_result.h"

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

    index_t parent = k;
    //index_t parent = parents_[k];
    index_t site = sites_[parent];
    avro_assert_msg( site < delaunay_.nb() , "elem = %lu, requested site %lu but nb_delaunay = %lu",k,site,delaunay_.nb());
    const real_t* z = delaunay_[ sites_[parent] ];
    T f = 0;
    for (coord_t d=0;d<delaunay_.dim()-1;d++)
      f += (z[d] - x[d])*(z[d] - x[d]);
    f -= weights_[ sites_[parent] ];
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
  volume_(sites.nb(),0.0),
  delaunay_(sites.dim()+1),
  exact_(true),
  diagram_(nullptr),
  decomposition_(decomposition_points_,domain.number()),
  decomposition_points_( points_.dim() ),
  facets_(points_,domain.number()-1)
{
  avro_assert( sites.nb() > 0 );
  coord_t dim = sites.dim();

  // copy the domain, considering that it is embedded into dim + 1
  std::vector<real_t> x( dim+1 , 0.0 );
  for (index_t k=0;k<domain.points().nb();k++)
  {
    for (index_t d=0;d<dim;d++)
      x[d] = domain.points()[k][d];
    domain_points_.create(x.data());
  }
  domain_.TopologyBase::copy( domain );

  // compute the maximum weight
  if (weight_.size() == 0)
    weight_.resize( sites.nb() , 0. );
  real_t wm = * std::max_element( weight_.begin() , weight_.end() ) + 1e-16;

  for (index_t k=0;k<sites.nb();k++)
  {
    for (coord_t d=0;d<delaunay_.dim();d++)
      x[d] = sites[k][d];
    delaunay_.create( x.data() );
    delaunay_[k][dim] = std::sqrt( wm - weight_[k] );
  }

  mass_.resize( delaunay_.nb() , 0.0 );
}

void
LaguerreDiagram::clear_decomposition()
{
  decomposition_points_.clear();
  decomposition_.clear();
}

class VoronoiFacets : public Topology<Polytope>
{
public:
  VoronoiFacets( Points& points , index_t number , const std::set<int>& bisectors , const std::map<SymbolicVertex,index_t>& symbolic ) :
    Topology<Polytope>(points,number),
    bisectors_(bisectors),
    symbolic_(symbolic)
  {
  }

  void extract()
  {
    std::vector< std::set<int> > bisector_vrep( bisectors_.size() );
    std::map<int,index_t> bisector_index;
    for (std::set<int>::const_iterator ib=bisectors_.begin();ib!=bisectors_.end();++ib)
      bisector_index.insert( {*ib,bisector_index.size()} );

    for (std::map<SymbolicVertex,index_t>::const_iterator iv=symbolic_.begin();iv!=symbolic_.end();++iv)
    {
      const SymbolicVertex& v = iv->first;
      const std::vector<int>& b = v.indices;

      for (index_t j=0;j<b.size();j++)
      {
        bisector_vrep[ bisector_index.at(b[j]) ].insert( iv->second );
      }
    }

    for (index_t k=0;k<bisector_vrep.size();k++)
    {
      std::vector<index_t> facet( bisector_vrep[k].begin(),bisector_vrep[k].end() );
      if (facet.size() == 0) continue;
      add( facet.data(),facet.size() );
    }

    //Table<index_t>::print();
    printf("found %lu facets\n",nb());
  }

private:
  const std::set<int>& bisectors_;
  const std::map<SymbolicVertex,index_t>& symbolic_;

  std::vector<int> elemL_;
  std::vector<int> elemR_;
};

void
LaguerreDiagram::compute()
{
  // compute the voronoi diagram
  if (diagram_ == nullptr)
    diagram_ = std::make_shared<VoronoiDiagram>( delaunay_ , domain_ , true );
  else
  {
    diagram_->clear();
    diagram_->points().clear();
  }
  diagram_->compute(exact_);

  time_neighbours_ = diagram_->time_neighbours();
  time_voronoi_    = diagram_->time_voronoi();

  coord_t dim = points_.dim();
  points_.clear();
  points_.set_dim(dim);
  clear();
  facets_.clear();

  // copy the points, removing the extra coordinate (in dim+1)
  for (index_t k=0;k<diagram_->points().nb();k++)
  {
    points_.create( diagram_->points()[k] );
    std::vector<int> bisectors = diagram_->points().incidence().get(k);
    points_.incidence().add( bisectors.data() , bisectors.size() );
  }
  TopologyBase::copy( *diagram_.get() );
  sites_ = diagram_->sites();

  avro_assert( sites_.size() == nb() );

  // calculate the set of all bisectors and symbolic vertices
  std::set<int> bisectors;
  std::map<SymbolicVertex,index_t> symbolic;
  const std::vector<SymbolicVertex>& vertices = diagram_->vertices();
  avro_assert( vertices.size() == points_.nb() );
  for (index_t k=0;k<points_.nb();k++)
  {
    std::vector<int> b = points_.incidence().get(k);
    for (index_t j=0;j<b.size();j++)
      bisectors.insert( b[j] );

    const SymbolicVertex& v = vertices[k];
    if (symbolic.find(v) == symbolic.end())
      symbolic.insert({v,k});
  }
  printf("there are %lu unique voronoi vertices (total %lu) and %lu bisectors (%lu)\n",symbolic.size(),points_.nb(),bisectors.size(),diagram_->bisectors().size());
  VoronoiFacets facets(points_,number_,bisectors,symbolic);
  facets.extract();
}

real_t
LaguerreDiagram::eval_objective( std::vector<real_t>& dE_dZ , std::vector<real_t>& dE_dW , std::vector<real_t>& V )
{
  // initialize size of centroids
  coord_t dim = points_.dim();
  Points centroids( dim );
  std::vector<real_t> x0( dim , 0. );
  for (index_t k=0;k<delaunay_.nb();k++)
    centroids.create( x0.data() );

  // size the volume array
  V.resize( delaunay_.nb() , 0. );

  real_t t0 = clock();
  #if 0
  // slow!
  SimplicialDecomposition<Polytope> decomposition(*this);
  decomposition.extract();
  real_t time_decompose = real_t(clock()-t0)/real_t(CLOCKS_PER_SEC);

  std::vector<index_t> simplices;
  std::vector<index_t> parents;
  decomposition.get_simplices( number() , simplices , parents );

  std::vector<index_t> sites( simplices.size()/(number_+1) );
  for (index_t k=0;k<sites.size();k++)
    sites[k] = sites_[parents[k]];

  #else
  // assumes cells have been decomposed into simplices
  std::vector<index_t> parents;
  std::vector<index_t> simplices;
  std::vector<index_t> sites;

  parents.reserve(1e6);
  simplices.reserve(1e6);
  sites.reserve(1e6);

  std::vector<index_t> simplex(number_+1);
  for (index_t k=0;k<delaunay_.nb();k++)
  {
    const VoronoiCell& cell_k = diagram_->cell(k);
    const Topology<Simplex>& cell_k_simplices = cell_k.simplices();

    index_t offset = decomposition_points_.nb();
    for (index_t j=0;j<cell_k_simplices.points().nb();j++)
      decomposition_points_.create( cell_k_simplices.points()[j] );

    for (index_t j=0;j<cell_k_simplices.nb();j++)
    {
      sites.push_back(k);
      parents.push_back(decomposition_.nb());
      for (coord_t d=0;d<cell_k_simplices.nv(j);d++)
        simplex[d] = cell_k_simplices(j,d) + offset;
      decomposition_.add( simplex.data() , simplex.size() );
    }
  }
  decomposition_.orient();
  avro_assert( fabs(decomposition_.volume() - 1.0) < 1e-12 ); // for now until more complicated geometries are studied
  index_t count = 0;
  simplices.resize( decomposition_.nb()*(number_+1) );
  for (index_t k=0;k<decomposition_.nb();k++)
  for (index_t j=0;j<decomposition_.nv(k);j++)
    simplices[count++] = decomposition_(k,j);
  real_t time_decompose = real_t(clock()-t0)/real_t(CLOCKS_PER_SEC);
  #endif
  index_t nb_simplices = simplices.size()/(number_+1);
  Topology<Simplex> simplex_topology( decomposition_.points() , number_ );
  const Simplex& element = decomposition_.element();
  for (index_t k=0;k<nb_simplices;k++)
  {
    simplex_topology.add( &simplices[k*(number_+1)] , number_+1 );

    real_t vk = element.volume( decomposition_.points() , &simplices[k*(number_+1)] , number_+1 );
    avro_assert( vk > 0. );

    index_t s = sites[k];
    avro_assert( s < V.size() );
    V[s] += vk;

    // compute the centroid of this piece
    std::vector<real_t> c( dim );
    numerics::centroid( &simplices[k*(number_+1)] , number_+1 , decomposition_.points() , c );

    // add the contribution to x*V and V
    for (coord_t d=0;d<centroids.dim();d++)
      centroids[s][d] += c[d]*vk;
  }

  // go back and divide by the volume
  for (index_t k=0;k<centroids.nb();k++)
  for (coord_t d=0;d<centroids.dim();d++)
  {
    // set the centroid to the original seed if the cell volume is zero
    if (V[k]<1e-12)
    {
      V[k] = 0.0;
      centroids[k][d] = delaunay_[k][d];
    }
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

    dE_dW[k] = mass_[k] -V[k];
    gnorm_w  += std::pow(dE_dW[k],2);
  }

  gnorm_x = std::sqrt( gnorm_x );
  gnorm_w = std::sqrt( gnorm_w );

  // compute the CVT energy
  ConicalProductQuadrature quadrature(number_,2);
  simplex_topology.element().set_basis( BasisFunctionCategory_Lagrange );
  quadrature.define();
  simplex_topology.element().load_quadrature(quadrature);

  typedef Integrand_Laguerre_Energy<real_t> Integrand_t;
  Integrand_t integrand(delaunay_,sites,parents,weight_,mass_);

  t0 = clock();
  Functional<Integrand_t> f(integrand);
  f.integrate( simplex_topology );
  real_t time_integrate = real_t( clock() - t0 )/real_t(CLOCKS_PER_SEC);

  real_t energy = f.value();
  for (index_t k=0;k<delaunay_.nb();k++)
    energy += weight_[k]*mass_[k];

  printf("--> f = %1.3e, |g_x| = %1.3e, |g_w| = %1.3e: t_n = %g sec. t_v = %g sec, t_d = %g sec., t_i = %g sec.\n",
          energy,gnorm_x,gnorm_w,time_neighbours_,time_voronoi_,time_decompose,time_integrate);

  return energy;
}

struct nlopt_data_cvt
{
	LaguerreDiagram& laguerre;
	index_t eval_count;
	real_t objective;
  int mode; // 0 for cvt, 1 for otm
};

real_t
LaguerreDiagram::eval_objective()
{
  std::vector<real_t> dE_dZ,dE_dW,volumes;
  return eval_objective(dE_dZ,dE_dW,volumes);
}

void
LaguerreDiagram::set_delaunay( const real_t* x , coord_t dim )
{
  for (index_t k=0;k<delaunay_.nb();k++)
  for (coord_t d=0;d<dim;d++)
    delaunay_[k][d] = x[k*dim+d];
}

void
LaguerreDiagram::set_weights( const real_t* w )
{
  for (index_t k=0;k<delaunay_.nb();k++)
    weight_[k] = w[k];
  coord_t dim = points_.dim();
  real_t wm = * std::max_element( weight_.begin() , weight_.end() );
  for (index_t k=0;k<weight_.size();k++)
    delaunay_[k][dim] = std::sqrt( wm - weight_[k] );
}

double
nlopt_otm_objective( unsigned n , const double* x , double* grad, void* data0 )
{
  nlopt_data_cvt* data = (nlopt_data_cvt*)(data0);

  LaguerreDiagram& laguerre = data->laguerre;

  printf("iteration %lu:\n",data->eval_count);
  data->eval_count++;

  // set the coordinates of the delaunay points
  if (data->mode == 0)
    laguerre.set_delaunay(x,laguerre.points().dim());
  else if (data->mode == 1)
    laguerre.set_weights(x);
  else
    avro_assert_not_reached;

  // recompute the voronoi diagram
  laguerre.compute();

  laguerre.clear_decomposition();

  // compute the energy and gradients
  std::vector<real_t> dE_dZ,dE_dW,volumes;
  real_t energy = laguerre.eval_objective( dE_dZ , dE_dW , volumes );
  if (grad != nullptr)
  {
    if (data->mode == 0)
    {
      // cvt mode
      for (index_t k=0;k<dE_dZ.size();k++)
        grad[k] = dE_dZ[k];
    }
    else if (data->mode == 1)
    {
      energy *= -1.0;

      // otm mode
      for (index_t k=0;k<dE_dW.size();k++)
        grad[k] = -dE_dW[k];
    }
    else
      avro_assert_not_reached;
  }
  laguerre.set_volumes(volumes);
  return energy;
}

#define OPTIMIZER_NLOPT 1
LaguerreDiagram* __laguerre__ = nullptr;
int __mode__ = -1;

void hlbfgs_otm_objective(int N, double* x, double *prev_x, double* f, double* g)
{
  avro_assert( __laguerre__ != nullptr );
  avro_assert( __mode__ >= 0 );
  LaguerreDiagram& laguerre = *__laguerre__;
  int mode = __mode__;

  if (mode == 0)
    laguerre.set_delaunay(x,laguerre.points().dim());
  else if (mode == 1)
    laguerre.set_weights(x);
  else
    avro_assert_not_reached;

  // recompute the voronoi diagram
  laguerre.compute();

  laguerre.clear_decomposition();

  // compute the energy and gradients
  std::vector<real_t> dE_dZ,dE_dW,volumes;
  real_t energy = laguerre.eval_objective( dE_dZ , dE_dW , volumes );
  if (g != nullptr)
  {
    if (mode == 0)
    {
      // cvt mode
      for (index_t k=0;k<dE_dZ.size();k++)
        g[k] = dE_dZ[k];
    }
    else if (mode == 1)
    {
      energy *= -1.0;

      // otm mode
      for (index_t k=0;k<dE_dW.size();k++)
        g[k] = -dE_dW[k];
    }
    else
      avro_assert_not_reached;
  }
  laguerre.set_volumes(volumes);

  *f = energy;
}

void newiteration(int iter, int call_iter, double *x, double* f, double *g,  double* gnorm)
{}

void
LaguerreDiagram::optimize_cvt()
{
  coord_t dim = points_.dim();
  index_t n = delaunay_.nb()*dim;

  std::vector<real_t> x(n);
  for (index_t k=0;k<delaunay_.nb();k++)
  for (index_t d=0;d<dim;d++)
    x[k*dim+d] = delaunay_(k,d);

  #if OPTIMIZER_NLOPT
  // setup the optimizer
  nlopt::opt opt( nlopt::LD_LBFGS , n );
  nlopt_data_cvt data = {*this,0,1,0};

	// set the objective function
	opt.set_min_objective( &nlopt_otm_objective , static_cast<void*>(&data) );

  // set some optimization parameters
  opt.set_xtol_rel(1e-12);
  opt.set_ftol_rel(1e-12);
  opt.set_maxeval(20);
  //opt.set_vector_storage(20);

  real_t f_opt;
  nlopt::result result = nlopt::result::SUCCESS;
  result = opt.optimize(x,f_opt);
  std::string desc = nloptResultDescription(result);
  printf("nlopt result: %s\n",desc.c_str());
  #else
  __laguerre__ = this;
  __mode__ = 0;
  double parameter[20];
  int info[20];
  int T = 0;
  int M = 7;
  int num_iter = 20;
  bool with_hessian = false;

  //initialize
  INIT_HLBFGS(parameter, info);
  info[4] = num_iter;
  info[6] = T;
  info[7] = with_hessian ? 1:0;
  info[10] = 0;
  info[11] = 1;

  HLBFGS(n, M, x.data() , hlbfgs_otm_objective , 0, HLBFGS_UPDATE_Hessian, newiteration, parameter, info);
  #endif
}

void
LaguerreDiagram::optimize_otm()
{
  index_t n = weight_.size();

  mass_.resize( n , real_t(1.0/weight_.size()) );
  for (index_t k=0;k<n;k++)
  {
    weight_[k] = 0;
  }

  mass_ = volume_;
  set_weights(weight_.data());

  real_t extra = 0;
  index_t f = 5;
  for (index_t k=0;k<10;k++)
  {
    extra += (f - 1)*mass_[k];
    mass_[k] *= f;
  }

  extra /= (mass_.size() - 10);
  for (index_t k=10;k<n;k++)
  {
    mass_[k] -= extra;
    //if (mass_[k] < 0) mass_[k] += extra;
  }

  std::vector<real_t> x(weight_.begin(),weight_.end());

  // setup the optimizer
  #if OPTIMIZER_NLOPT
  nlopt::opt opt( nlopt::LD_LBFGS , n );
  nlopt_data_cvt data = {*this,0,1,1};

	// set the objective function
	opt.set_min_objective( &nlopt_otm_objective , static_cast<void*>(&data) );

  // set some optimization parameters
  opt.set_xtol_rel(1e-12);
  opt.set_ftol_rel(1e-12);
  opt.set_maxeval(20);
  //opt.set_vector_storage(20);

  // set the lower and upper bounds on the entries of the step matrix
	std::vector<real_t> lower_bound( n , 0.0 );
	std::vector<real_t> upper_bound( n ,  1e10 );
	opt.set_lower_bounds(lower_bound);
	opt.set_upper_bounds(upper_bound);

  real_t f_opt;
  nlopt::result result = nlopt::result::SUCCESS;
  result = opt.optimize(x,f_opt);
  std::string desc = nloptResultDescription(result);
  printf("nlopt result: %s\n",desc.c_str());
  #else

  __laguerre__ = this;
  __mode__ = 1;
  double parameter[20];
  int info[20];
  int T = 0;
  int M = 7;
  int num_iter = 20;
  bool with_hessian = false;

  //initialize
  INIT_HLBFGS(parameter, info);
  info[4] = num_iter;
  info[6] = T;
  info[7] = with_hessian ? 1:0;
  info[10] = 0;
  info[11] = 1;

  HLBFGS(n, M, x.data() , hlbfgs_otm_objective , 0, HLBFGS_UPDATE_Hessian, newiteration, parameter, info);
  #endif
}

} // delaunay

} // avro
