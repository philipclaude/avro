#include "common/error.h"
#include "common/process.h"

#include "element/quadrature.h"

#include "mesh/decomposition.h"

#include "numerics/geometry.h"
#include "numerics/integration.h"
#include "numerics/nlopt_result.h"

#include "voronoi/power.h"
#include "voronoi/voronoi_cell.h"

#include <HLBFGS/HLBFGS.h>
#include <HLBFGS/Lite_Sparse_Matrix.h>

#include <tinymat/types/SurrealS.h>
#include <tinymat/types/SurrealD.h>

#include <nlopt.hpp>

#include <iomanip>

namespace avro
{

namespace delaunay
{

template<typename _T>
class Integrand_Power_Energy : public Integrand<Integrand_Power_Energy<_T>>
{
public:

  typedef _T T;

  Integrand_Power_Energy( const Points& delaunay , const std::vector<index_t>& sites , const std::vector<index_t>& parents ,
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
    index_t site = sites_[parent];
    avro_assert_msg( site < delaunay_.nb() , "elem = %lu, requested site %lu but nb_delaunay = %lu",k,site,delaunay_.nb());
    const real_t* z = delaunay_[ sites_[parent] ];
    T f = 0;
    for (coord_t d=0;d<delaunay_.dim()-1;d++)
      f += (z[d] - x[d])*(z[d] - x[d]);
    return f;
  }

  private:
    const Points& delaunay_;
    const std::vector<index_t>& sites_;
    const std::vector<index_t>& parents_;
    const std::vector<T>& weights_;
    const std::vector<real_t>& mass_;
};

PowerDiagram::PowerDiagram( const Points& sites , const Topology<Simplex>& domain , const std::vector<real_t>& weights ) :
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
  facets_(*this)
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
PowerDiagram::clear_decomposition()
{
  decomposition_points_.clear();
  decomposition_.clear();
}

PowerFacets::PowerFacets( PowerDiagram& diagram ) :
  Topology<Polytope>(diagram.points(),diagram.number()-1),
  diagram_(diagram),
  bij_( diagram.points().dim() )
{}

void
PowerFacets::extract()
{
  // clear any existing information
  bisectors_.clear();
  symbolic_.clear();
  unique_vertex_.clear();
  label_.clear();
  cellR_.clear();
  cellL_.clear();
  clear();

  // calculate the set of all bisectors and symbolic vertices
  const std::vector<SymbolicVertex>& vertices = diagram_.diagram().symbolic_vertices();
  avro_assert( vertices.size() == points_.nb() );
  for (index_t k=0;k<points_.nb();k++)
  {
    // retrieve all bisectors this point is on
    std::vector<int> b = points_.incidence().get(k);
    for (index_t j=0;j<b.size();j++)
      bisectors_.insert( b[j] );

    // associate this ver
    const SymbolicVertex& v = vertices[k];
    if (symbolic_.find(v) == symbolic_.end())
    {
      symbolic_.insert({v,k});
      unique_vertex_.insert( {k,k} );
    }
    else
      unique_vertex_.insert( {k,symbolic_.at(v) } );
  }
  /*printf("there are %lu unique voronoi vertices (total %lu) and %lu bisectors (%lu)\n",
          symbolic_.size(),points_.nb(),bisectors_.size(),diagram_.diagram().bisectors().size());*/

  std::vector<int> hrep;

  std::map<ElementIndices,index_t> created;
  ElementIndices facet;
  for (index_t k=0;k<diagram_.nb();k++)
  {
    // extract the hrep of this cell
    hrep.clear();
    diagram_.element().hrep( diagram_(k) , diagram_.nv(k) , hrep );

    // loop through the hrep
    for (index_t j=0;j<hrep.size();j++)
    {
      // extract the vrep and determine unique entries
      std::vector<index_t>& vrep = facet.indices;
      vrep.clear();
      diagram_.element().vrep( diagram_(k) , diagram_.nv(k) , hrep[j] , vrep );
      for (index_t i=0;i<vrep.size();i++)
        vrep[i] = unique_vertex_.at( vrep[i] );
      std::sort( vrep.begin(),vrep.end());

      // determine if this facet was already created
      if (created.find(facet) != created.end())
      {
        // created: assign the right cell
        cellR_[created.at(facet)] = k;
        avro_assert( label_[created.at(facet)] == hrep[j] );
      }
      else
      {
        // not created: assign the left cell information
        cellL_.push_back(k);
        cellR_.push_back(-1);
        label_.push_back(hrep[j]);
        created.insert({facet,nb()});
        add( vrep.data() , vrep.size() );
      }
    }
  }
  compute_quantities();
}

void
PowerFacets::compute_quantities()
{
  bij_.clear();
  aij_.clear();
  neighbours_.clear();
  neighbours_.resize( diagram_.delaunay().nb() );
  indices_.clear();

  // first determine how many quantities we need
  std::map<index_t,Bisector> facet2bisector;
  std::map<Bisector,index_t> bisector2facet;
  std::set<Bisector> bisectors;
  for (index_t k=0;k<nb();k++)
  {
    // skip interior domain facets
    if (label_[k] < 0 && cellR_[k] >= 0) continue;
    if (label_[k] < 0)
    {
      // this is a boundary facet
      // we need to determine the site of the left element
      avro_assert( cellR_[k] < 0 );

      index_t pi = diagram_.site( cellL_[k] );
      Bisector b(label_[k],pi);
      if (bisectors.find(b) == bisectors.end()) bisectors.insert(b);
      facet2bisector.insert({k,b});
      bisector2facet.insert({b,k});
      neighbours_[pi].insert(label_[k]);
      continue;
    }
    const Bisector& b = diagram_.diagram().bisectors().at( label_[k] );
    if (bisectors.find(b) == bisectors.end()) bisectors.insert(b);
    facet2bisector.insert({k,b});
    bisector2facet.insert({b,k});
    neighbours_[b.p0].insert(b.p1);
    neighbours_[b.p1].insert(b.p0);
  }
  avro_assert( bisector2facet.size() == bisectors.size() );

  std::map<index_t,index_t> facet2id;
  std::map<index_t,index_t> id2facet;

  aij_.resize( bisectors.size() );
  std::vector<real_t> p0(points_.dim(),0.0);
  std::map<Bisector,index_t>::iterator it;
  for (it=bisector2facet.begin();it!=bisector2facet.end();++it)
  {
    index_t id = bij_.nb();
    facet2id.insert( {it->second,id} );
    id2facet.insert( {id,it->second} );
    bij_.add( p0.data() , p0.size() );
    aij_[id] = 0.;
    indices_.insert( {id,it->first} );
  }

  std::vector<real_t> centroid(points_.dim());
  for (it=bisector2facet.begin();it!=bisector2facet.end();++it)
  {
    index_t facet = it->second;
    index_t id = facet2id[facet];

    real_t volume = 0.0;
    if (number_ == 1)
    {
      // compute the volume of this facet
      avro_assert( nv(facet) == 2 );
      volume = numerics::distance( points_[operator()(facet,0)] , points_[operator()(facet,1)] , points_.dim() );
      aij_[id] += volume;
    }
    else
      avro_implement;

    // compute the centroid of this facet
    numerics::centroid( operator()(facet) , nv(facet) , points_ , centroid );

    // add the contribution (we will divide by the total volume later)
    for (index_t j=0;j<centroid.size();j++)
      bij_[id][j] += volume*centroid[j];
  }

  for (index_t k=0;k<bij_.nb();k++)
  for (index_t j=0;j<points_.dim();j++)
    bij_[k][j] /= aij_[k];
}

void
PowerFacets::print() const
{
  for (index_t k=0;k<nb();k++)
  {
    std::vector<index_t> f = get(k);
    print_inline( f , "cellL = " + std::to_string(cellL_[k]) + ", cellR = " + std::to_string(cellR_[k]) + ", bisector = " + std::to_string(label_[k]) + " facet =" );
  }
}

class Particles : public Topology<Polytope>
{
public:
  Particles( const PowerDiagram& power_diagram )
  {
    extract(power_diagram);
  }

  // extraction
  void extract(const PowerDiagram& power_diagram , bool unique_vertices=true )
  {
    clear();

    // construct the list of bisectors and symbolic vertices

    // create the vertices (option to be unique based on the symbolic info)

    // add the cells, carefully treating vertices on mesh facets

  }

  void clear()
  {
    Topology<Polytope>::clear();
    points_.clear();
    symbolic_.clear();
    bisectors_.clear();
  }

private:
  Points points_;
  std::vector<SymbolicVertex> symbolic_;
  std::map<int,Bisector> bisectors_;
};

void
PowerDiagram::compute()
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

  facets_.extract();
}


index_t __iter__ = 0;
int __mode__ = -1;

real_t
PowerDiagram::eval_objective( std::vector<real_t>& dE_dZ , std::vector<real_t>& dE_dW , std::vector<real_t>& V )
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

    dE_dW[k] = -(V[k] - mass_[k]);
    gnorm_w  += std::pow(dE_dW[k],2);
  }

  gnorm_x = std::sqrt( gnorm_x );
  gnorm_w = std::sqrt( gnorm_w );

  // compute the CVT energy
  ConicalProductQuadrature quadrature(number_,2);
  simplex_topology.element().set_basis( BasisFunctionCategory_Lagrange );
  quadrature.define();
  simplex_topology.element().load_quadrature(quadrature);

  typedef Integrand_Power_Energy<real_t> Integrand_t;
  Integrand_t integrand(delaunay_,sites,parents,weight_,mass_);

  t0 = clock();
  Functional<Integrand_t> f(integrand);
  f.integrate( simplex_topology );
  real_t time_integrate = real_t( clock() - t0 )/real_t(CLOCKS_PER_SEC);

  real_t energy = f.value();
  for (index_t k=0;k<delaunay_.nb();k++)
    energy -= weight_[k]*(V[k] - mass_[k]);

  //printf("--> f = %1.3e, |g_x| = %1.3e, |g_w| = %1.3e: t_n = %g sec. t_v = %g sec\n",
  //        energy,gnorm_x,gnorm_w,time_neighbours_,time_voronoi_);

  UNUSED(time_integrate);
  UNUSED(time_decompose);

  real_t gnorm = (__mode__ == 0) ? gnorm_x : gnorm_w;

  std::cout <<
  std::setw(6) << std::left << __iter__++ << "|" <<
  std::setprecision(4) << std::left << std::scientific << std::fabs(energy) << "|" <<
  std::setprecision(4) << std::left << std::scientific << gnorm << "|" <<
  std::setw(11) << std::right << std::fixed << time_voronoi_ << "|" <<
  std::setw(11) << std::right << std::fixed << time_neighbours_ <<
  std::endl;

  return energy;
}

void
PowerDiagram::get_hessian( Lite_Sparse_Matrix& h ) const
{
  h.begin_fill_entry();

  // loop through all the facets
  const std::vector<real_t> aij = facets_.aij();
  const DOF<real_t>& bij = facets_.bij();
  const std::map<index_t,Bisector>& bisectors = facets_.indices();

  avro_assert( aij.size() == bij.nb() );
  avro_assert( aij.size() == bisectors.size() );

  std::vector<real_t> d( delaunay_.nb() , 0.0 );

  for (index_t k=0;k<aij.size();k++)
  {
    const Bisector& b = bisectors.at(k);
    avro_assert( b.p0 >= 0 );


    if (b.p1 < 0)
    {
      // skip boundary facets (we apply neumann boundary conditions)
      continue;
    }

    index_t r = b.p0;
    index_t c = b.p1;

    real_t lij = numerics::distance( delaunay_[r] , delaunay_[c] , delaunay_.dim() );
    real_t eij = -0.5*aij[k]/lij;

    h.fill_entry(r,c, eij );
    h.fill_entry(c,r, eij );
    d[r] += eij;
    d[c] += eij;
  }

  for (index_t k=0;k<delaunay_.nb();k++)
    h.fill_entry(k,k , d[k] );

  h.end_fill_entry();
}

struct nlopt_data_cvt
{
	PowerDiagram& power;
	index_t eval_count;
	real_t objective;
  int mode; // 0 for cvt, 1 for otm
};

void
PowerDiagram::set_delaunay( const real_t* x , coord_t dim )
{
  for (index_t k=0;k<delaunay_.nb();k++)
  for (coord_t d=0;d<dim;d++)
    delaunay_[k][d] = x[k*dim+d];
}

void
PowerDiagram::set_weights( const real_t* w )
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

  PowerDiagram& power = data->power;

  //printf("iteration %lu:\n",data->eval_count);
  data->eval_count++;

  // set the coordinates of the delaunay points
  if (data->mode == 0)
    power.set_delaunay(x,power.points().dim());
  else if (data->mode == 1)
    power.set_weights(x);
  else
    avro_assert_not_reached;

  // recompute the voronoi diagram
  power.compute();
  power.clear_decomposition();

  // compute the energy and gradients
  std::vector<real_t> dE_dZ,dE_dW,volumes;
  real_t energy = power.eval_objective( dE_dZ , dE_dW , volumes );
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
  power.set_volumes(volumes);
  return energy;
}

#define OPTIMIZER_NLOPT 1
PowerDiagram* __power__ = nullptr;
Lite_Sparse_Matrix* __h_matrix__ = nullptr;

void hlbfgs_otm_objective(int N, double* x, double *prev_x, double* f, double* g)
{
  avro_assert( __power__ != nullptr );
  avro_assert( __mode__ >= 0 );
  PowerDiagram& power = *__power__;
  int mode = __mode__;

  if (mode == 0)
    power.set_delaunay(x,power.points().dim());
  else if (mode == 1)
    power.set_weights(x);
  else
    avro_assert_not_reached;

  // recompute the voronoi diagram
  power.compute();
  power.clear_decomposition();

  // compute the energy and gradients
  std::vector<real_t> dE_dZ,dE_dW,volumes;
  real_t energy = power.eval_objective( dE_dZ , dE_dW , volumes );
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
  power.set_volumes(volumes);

  *f = energy;
}

void hlbfgs_otm_objective_hessian(int N, double* x, double *prev_x, double* f, double* g, HESSIAN_MATRIX& h)
{
  avro_assert( __power__ != nullptr );
  avro_assert( __mode__ >= 0 );
  PowerDiagram& power = *__power__;
  int mode = __mode__;

  if (mode == 0)
    power.set_delaunay(x,power.points().dim());
  else if (mode == 1)
    power.set_weights(x);
  else
    avro_assert_not_reached;

  // recompute the voronoi diagram
  power.compute();
  power.clear_decomposition();

  // compute the energy and gradients
  std::vector<real_t> dE_dZ,dE_dW,volumes;
  real_t energy = power.eval_objective( dE_dZ , dE_dW , volumes );
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
      // otm mode: energy is concave so make this a minimization problem
      energy *= -1.0;
      for (index_t k=0;k<dE_dW.size();k++)
        g[k] = -dE_dW[k];

      // add the hessian
      if (__h_matrix__ != nullptr) delete __h_matrix__;
      __h_matrix__ = new Lite_Sparse_Matrix(N,N,SYM_BOTH,CCS,FORTRAN_TYPE,true);
      power.get_hessian(*__h_matrix__);

      h.set_diag(__h_matrix__->get_diag());
      h.set_values(__h_matrix__->get_values());
      h.set_rowind(__h_matrix__->get_rowind());
      h.set_colptr(__h_matrix__->get_colptr());
      h.set_nonzeros(__h_matrix__->get_nonzero());
    }
    else
      avro_assert_not_reached;
  }
  power.set_volumes(volumes);

  *f = energy;
}

void newiteration(int iter, int call_iter, double *x, double* f, double *g,  double* gnorm)
{}

void
PowerDiagram::optimize_cvt()
{
  __iter__ = 0;
  printf("-----------------------------------------------------\n");
  printf("%5s|%10s|%10s|%10s|%10s\n","iter  "," energy "," gradient "," t_vor (s) "," t_knn (s) ");
  printf("-----------------------------------------------------\n");

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
  __power__ = this;
  __mode__ = 0;
  double parameter[20];
  int info[20];
  int T = 1;
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
PowerDiagram::optimize_otm()
{
  __iter__ = 0;
  printf("-----------------------------------------------------\n");
  printf("%5s|%10s|%10s|%10s|%10s\n","iter  "," energy "," gradient "," t_vor (s) "," t_knn (s) ");
  printf("-----------------------------------------------------\n");

  index_t n = weight_.size();

  mass_.resize( n , real_t(1.0/weight_.size()) );
  for (index_t k=0;k<n;k++)
  {
    weight_[k] = 0;
  }

  mass_ = volume_;
  set_weights(weight_.data());

  for (index_t k=0;k<n;k++)
  {
    for (coord_t d=0;d<delaunay_.dim();d++)
    {
      real_t vd = random_within(-0.01,0.01);
      real_t xd = delaunay_[k][d] + vd;
      if (xd < 0 || xd > 1) continue;
      delaunay_[k][d] += vd;
    }
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

  __power__ = this;
  __mode__ = 1;
  double parameter[20];
  int info[20];
  int T = 0;
  int M = 7;
  int num_iter = 20;
  bool with_hessian = true;

  // initialize
  INIT_HLBFGS(parameter, info);
  parameter[5] = 1e-12;
  info[3] = 1;
  info[4] = num_iter;
  info[6] = T;
  info[5] = 1; // verbose = 1
  info[7] = with_hessian ? 1:0;
  info[10] = 0;
  info[11] = 1;

  if (with_hessian)
    HLBFGS(n, M, x.data() , hlbfgs_otm_objective , hlbfgs_otm_objective_hessian, HLBFGS_UPDATE_Hessian, newiteration, parameter, info);
  else
    HLBFGS(n, M, x.data() , hlbfgs_otm_objective , 0 , HLBFGS_UPDATE_Hessian, newiteration, parameter, info);
  if (__h_matrix__ != nullptr) delete __h_matrix__;
  #endif
}

/*
------------------------------------------------------------------------------------
HLBFGS reference: https://xueyuhanlang.github.io/software/HLBFGS/
------------------------------------------------------------------------------------
Entry	        Default value	Meaning
------------------------------------------------------------------------------------
PARAMETER[0]	1.0e-4	      function tolerance used in line-search
PARAMETER[1]	1.0e-16	      variable tolerance used in line-search
PARAMETER[2]	0.9	          gradient tolerance used in line-search
PARAMETER[3]	1.0e-20	      stpmin used in line-search
PARAMETER[4]	1.0e+20	      stpmax used in line-search
PARAMETER[5]	1.0e-9	      the stop criterion ( ||G||/max(1,||X||) < PARAMETER[5] )
PARAMETER[6]	1.0e-10	      the stop criterion ( ||G|| < PARAMETER[6] )
------------------------------------------------------------------------------------
INFO[0]	      20	          the max number of evaluation in line-search
INFO[1]	      0	            the total number of evalfunc calls
INFO[2]	      0	            the current number of iterations
INFO[3]	      0	            the lbfgs strategy. 0: standard, 1: M1QN3 strategy[8](recommended).
INFO[4]	      100000	      the max number of iterations
INFO[5]	      1	            1: print message, 0: do nothing
INFO[6]	      10	          T: the update interval of Hessian. (typical choices: 0-200)
INFO[7]	      0	            0: without hessian, 1: with accurate hessian
INFO[8]	      15	          icfs parameter
INFO[9]	      0	            0: classical line-search; 1: modified line-search (it is not useful in practice)
INFO[10]	    0	            0: Disable preconditioned CG; 1: Enable preconditioned CG
INFO[11]	    1	            0 or 1 defines different methods for choosing beta in CG.
INFO[12]	    1	            internal usage. 0: only update the diag in USER_DEFINED_HLBFGS_UPDATE_H; 1: default.
------------------------------------------------------------------------------------
In general the users only need to specify the values of INFO[3,4,5,6,7,10]. For different optimization methods, please check the following table.
------------------------------------------------------------------------------------
Method	                                          Setting
------------------------------------------------------------------------------------
Gradient Decent	                                  M=0, INFO[7]=0, INFO[10]=0
Conjugate Gradient	                              M=0, INFO[7]=0, INFO[10]=1
Newton's method	                                  M=0, INFO[6]=0, INFO[7]=1, INFO[10]=0
L-BFGS	                                          M>=1, INFO[3]=0, INFO[7]=0, INFO[10]=0
M1QN3	                                            M>=1, INFO[3]=1, INFO[7]=0, INFO[10]=0
Preconditioned L-BFGS	                            M>=1, INFO[3]=0, INFO[6]=T>=0, INFO[7]=1, INFO[10]=0
Preconditioned M1QN3                             	M>=1, INFO[3]=1, INFO[6]=T>=0, INFO[7]=1, INFO[10]=0
Preconditioned Conjugate Gradient without Hessian	M>=1, INFO[3]=0 or 1, INFO[7]=0, INFO[10]=1
Preconditioned Conjugate Gradient with Hessian	  M>=1, INFO[3]=0 or 1, INFO[6]=T>=0, INFO[7]=1, INFO[10]=1
------------------------------------------------------------------------------------
*/

} // delaunay

} // avro
