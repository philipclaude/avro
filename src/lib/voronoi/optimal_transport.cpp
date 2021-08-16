#include "common/parallel_for.h"

#include "element/quadrature.h"

#include "numerics/geometry.h"
#include "numerics/integration_rank.h"
#include "numerics/nlopt_result.h"

#include "voronoi/optimal_transport.h"

#include <nlopt.hpp>
#include "json/json.hpp"

#include <fstream>
#include <iomanip>

// used in the side5 predicate
bool __check_capacity__ = true;

namespace avro
{

namespace voronoi
{

template<typename type>
LaguerreDiagram<type>::LaguerreDiagram( Points& delaunay , const Topology<type>& domain ) :
  Topology<Polytope>(points_,domain.number()),
  points_(delaunay.dim()),
  delaunay_(delaunay),
  domain_(domain),
  neighbour_search_(nullptr)
{}

template<>
void
LaguerreDiagram<Simplex>::initialize()
{
  // we need the facets but not the edges
  domain_facets_ = std::make_shared<RVDFacets>(domain_);
}

template<>
void
LaguerreDiagram<Simplex>::create( bool exact , index_t nb_nns )
{
  // create a cell for every simplex in the mesh
  cells_.resize( domain_.nb() , nullptr );
  for (index_t k=0;k<domain_.nb();k++)
  {
    cells_[k] = std::make_shared<LaguerreCell<Simplex>>(k,delaunay_,*neighbour_search_,domain_,domain_facets_.get(),exact,nb_nns);
    cells_[k]->set_facets( domain_facets_.get() );
    cells_[k]->set_edges( domain_edges_ );
  }
}

template<>
void
LaguerreDiagram<Polytope>::create( bool exact , index_t nb_nns )
{
  // create a cell for every delaunay vertex
  cells_.resize( delaunay_.nb() , nullptr );
  for (index_t k=0;k<delaunay_.nb();k++)
  {
    if (cells_[k] == nullptr || neighbour_search_ == nullptr)
      cells_[k] = std::make_shared<LaguerreCell<Polytope>>(k,delaunay_,*neighbour_search_,domain_,exact,nb_nns);
    cells_[k]->set_facets( domain_facets_.get() );
    cells_[k]->set_edges( domain_edges_ );
  }
}

template<>
void
LaguerreDiagram<Polytope>::initialize()
{
  // we need the edges but not the facets
  domain_edges_.clear();
  domain_.get_edges(domain_edges_);
  elem_.resize( delaunay_.nb() , 0 );
}

template<typename type>
void
LaguerreDiagram<type>::compute( bool exact , IntegrationSimplices* triangulation )
{
  // initialize the domain data
  initialize();

  // initialize the nearest neighbours
  real_t t0 = clock();
  index_t nb_nns = 50;
  if (delaunay_.nb() < nb_nns) nb_nns = delaunay_.nb();
  const coord_t dim = delaunay_.dim();
  if (neighbour_search_ == nullptr)
    neighbour_search_ = GEO::NearestNeighborSearch::create(dim,"ANN");
  std::vector<real_t> x(delaunay_.nb()*dim);
  for (index_t k=0;k<delaunay_.nb();k++)
  for (index_t d=0;d<dim;d++)
    x[k*dim+d] = delaunay_[k][d];
  t0 = clock();
  neighbour_search_->set_points( delaunay_.nb() , x.data() );
  time_neighbours_ = real_t( clock() - t0 )/real_t(CLOCKS_PER_SEC);

  // create the laguerre cells
  create(exact,nb_nns);

  // clip the cells
  #if 1
  ProcessCPU::parallel_for(
    parallel_for_member_callback( this , &thisclass::clip ),
    0,cells_.size() );
  #else
  for (index_t k=0;k<cells_.size();k++) {
    cells_[k]->compute();
  }
  #endif

  time_voronoi_   = 0;
  time_decompose_ = 0;

  // accumulate the result
  sites_.clear();
  std::map<Bisector,int> bisectors;
  for (index_t k=0;k<cells_.size();k++) {

    // retrieve the voronoi cell
    const LaguerreCell<type>& cell = *cells_[k].get();

    // accumulate the timing information
    time_neighbours_ += cell.time_neighbours() / ProcessCPU::maximum_concurrent_threads();
    time_voronoi_    += cell.time_clip() / ProcessCPU::maximum_concurrent_threads();
    time_decompose_  += cell.time_decompose() / ProcessCPU::maximum_concurrent_threads();

    // cells could be submerged for weighted sites
    if (cell.nb() == 0) continue;

    #if 1 // turn off for lower memory usage

    // add the cells
    for (index_t j=0;j<cell.nb();j++) {
      sites_.push_back( k );
      std::vector<index_t> polytope = cell.get(j);

      for (index_t i=0;i<polytope.size();i++)
        polytope[i] += points_.nb();
      add( polytope.data() , polytope.size() );
    }

    // add the points
    for (index_t j=0;j<cell.points().nb();j++) {
      points_.create( cell.points()[j] );
    }

    // add the vertex-to-facet information
    const Table<int>& vf = cell.element().incidence();
    for (index_t j=0;j<vf.nb();j++) {

      std::vector<int> f = vf.get(j);

      // map the bisctors to a global numbering
      for (index_t j = 0; j < f.size(); j++) {
        if (f[j] < 0) continue; // mesh facets are in global numbering
        index_t pi, pj;
        cell.get_bisector( f[j] , pi , pj );
        Bisector b(pi,pj);
        int h;
        if (bisectors.find(b) == bisectors.end()) {
          h = bisectors.size();
          bisectors.insert({b,h});
        }
        else {
          h = bisectors.at(b);
        }
        f[j] = h;
      }

      points_.incidence().add( f.data() , f.size() );
    }
    #endif

    // option to add the triangulation data
    if (triangulation == nullptr) continue;

    const IntegrationSimplices& tk = cell.simplices();
    std::vector<index_t> simplex(tk.number()+1);
    for (index_t j = 0; j < tk.nb(); j++)
    {
      for (index_t i = 0; i < tk.nv(j); i++)
        simplex[i] = tk(j,i) + triangulation->points().nb();
      triangulation->add_simplex( simplex , tk.simplex2elem(j) , tk.simplex2site(j) );
    }

    // this will create points that have the same dimension as the requested triangulation
    // in other words, even though the laguerre cell might be in dim+1, we only copy dim coordinates
    for (index_t j = 0; j < tk.points().nb(); j++)
      triangulation->add_point( tk.points()[j] , tk.point2elem(j) , tk.point2site(j) );
  }
  //printf("--> timing: neighbours (%3.4g sec.), clipping (%3.4g sec.), decompose (%3.4g sec.)\n",time_neighbours_,time_voronoi_,time_decompose);
  compute_neighbour_properties();
}

template<typename type>
void
LaguerreDiagram<type>::compute_neighbour_properties()
{
  neighbour_counts_.clear();
  neighbour_counts_.resize( delaunay_.nb() , 0 );

  std::vector<int> hrep;
  for (index_t k = 0; k < cells_.size(); k++)
  {
    for (index_t j = 0; j < cells_[k]->nb(); j++)
    {
      hrep.clear();
      cells_[k]->element().hrep( (*cells_[k])(j) , cells_[k]->nv(j) , hrep );

      index_t count = 0;
      index_t site = cells_[k]->cell2site(j);
      for (index_t i = 0; i < hrep.size(); i++)
      {
        //if (hrep[i] < 0) continue; // mesh facet

        // decode the bisector associated with this?
        count++;
      }
      neighbour_counts_[site] += count;
    }
  }

  index_t total_neighbours = 0;
  for (index_t k = 0; k < cells_.size(); k++)
    total_neighbours += neighbour_counts_[k];

  average_neighbours_ = real_t(total_neighbours) / real_t(cells_.size());
  minimum_neighbours_ = * std::min_element( neighbour_counts_.begin() , neighbour_counts_.end() );
  maximum_neighbours_ = * std::max_element( neighbour_counts_.begin() , neighbour_counts_.end() );
}

template<typename type>
LaguerreDiagram<type>::~LaguerreDiagram() {
  if (neighbour_search_ != nullptr) delete neighbour_search_;
}

real_t
Integrand_Transport_Energy::operator()( index_t k , const QuadraturePoint& point , const real_t* x ) const
{
  avro_assert_msg( k < simplices_.nb() , "elem = %lu, nb_simplices = %lu", k , simplices_.nb() );

  // retrieve the delaunay site and coordinate
  index_t site = simplices_.simplex2site(k);
  avro_assert_msg( site < delaunay_.nb() , "elem = %lu, requested site %lu but nb_delaunay = %lu",k,site,delaunay_.nb());
  const real_t* z = delaunay_[site];

  // evaluate the density
  real_t rho = density_.evaluate( k , point.coordinate() , x );

  T f = 0;
  for (coord_t d=0;d<dim_;d++)
    f += rho*(z[d] - x[d])*(z[d] - x[d]);
  return f;
}

class Integrand_Mass : public Integrand<Integrand_Mass>
{
public:

  typedef real_t T;

  Integrand_Mass(const DensityMeasure& density) :
    density_(density)
  {}

  bool needs_solution() const { return false; }
  bool needs_gradient() const { return false; }
  bool needs_hessian() const { return false; }

  T operator()( index_t k , const QuadraturePoint& point , const real_t* x ) const
  {
    return density_.evaluate(k,point.coordinate(),x);
  }

private:
  const DensityMeasure& density_;
};

class Integrand_Moment : public Integrand<Integrand_Moment>
{
public:

  typedef real_t T;

  Integrand_Moment(coord_t dim, const DensityMeasure& density) :
    dim_(dim),
    density_(density)
  {}

  bool needs_solution() const { return false; }
  bool needs_gradient() const { return false; }
  bool needs_hessian() const { return false; }

  T operator()( index_t k , const QuadraturePoint& point , const real_t* x , std::vector<T>& I ) const
  {
    real_t rho = density_.evaluate(k,point.coordinate(),x);
    for (index_t r = 0; r < dim_; r++)
      I[r] = rho*x[r];
    return 1.0;
  }

private:
  coord_t dim_;
  const DensityMeasure& density_;
};

template<typename type>
SemiDiscreteOptimalTransport<type>::SemiDiscreteOptimalTransport( const Topology<type>& domain , DensityMeasure* density ) :
  domain_(domain),
  density_(density),
  delaunay_(domain.points().dim()),
  diagram_(delaunay_,domain_),
  simplices_(domain.number(),domain.number()),
  exact_(false),
  print_(true),
  quad_order_(4),
  weight_max_(1e6),
  save_every_(1e6),
  prefix_("avro-sdot")
{}

template<typename type>
void
SemiDiscreteOptimalTransport<type>::sample( index_t nb_samples )
{
  coord_t dim = domain_.points().dim();

  // create samples in a dim+1 space, but leave the last dimension as 0
  delaunay_.clear();
  for (index_t k=0;k<nb_samples;k++)
  {
    std::vector<real_t> p(dim,0.);
    for (coord_t d = 0; d < domain_.number(); d++)
      p[d] = random_within(0.0,1.0);
    delaunay_.create(p.data());
  }
  nu_.resize( delaunay_.nb() , 0.0 );
  weight_.resize( delaunay_.nb() , 0.0 );
}

template<typename type>
void
SemiDiscreteOptimalTransport<type>::compute_laguerre()
{
  // compute the laguerre diagram along with the integration simplices
  diagram_.clear();
  diagram_.points().clear();
  simplices_.clear();
  diagram_.compute(exact_,&simplices_);
}

static real_t
calculate_norm( const real_t* x , index_t nb )
{
  real_t n = 0.0;
  for (index_t k = 0; k < nb; k++)
    n += x[k]*x[k];
  return std::sqrt(n);
}

template<typename type>
real_t
SemiDiscreteOptimalTransport<type>::evaluate( real_t* dc_dx , real_t* dc_dw )
{
  avro_assert( density_ != nullptr );
  clock_t t0;
  real_t time_integrate = 0.0;

  coord_t dim = domain_.number();

  // compute the laguerre diagram
  compute_laguerre();

  // define the integration rules
  simplices_.element().set_basis( BasisFunctionCategory_Lagrange );

  // get the masses
  t0 = clock();
  typedef voronoi::Integrand_Mass Integrand_Mass_t;
  Integrand_Mass_t integrandm(*density_);
  Functional<Integrand_Mass_t> fm(integrandm,simplices_);
  std::vector<real_t> masses(simplices_.nb());
  fm.integrate(simplices_,masses.data());
  time_integrate += real_t( clock() - t0 )/real_t(CLOCKS_PER_SEC);

  // get the moments
  t0 = clock();
  typedef voronoi::Integrand_Moment Integrand_Moment_t;
  Integrand_Moment_t integrandxm(dim,*density_);
  Functional_Ranked<Integrand_Moment_t> fxm(integrandxm,dim,simplices_);
  std::vector<real_t> moments( dim * simplices_.nb() );
  fxm.integrate(simplices_,moments.data());
  time_integrate += real_t( clock() - t0 )/real_t(CLOCKS_PER_SEC);

  // first compute the mass of each voronoi cell
  mass_.clear();
  mass_.resize( delaunay_.nb() , 0.0 );
  for (index_t k = 0; k < simplices_.nb(); k++)
  {
    index_t s = simplices_.simplex2site(k);
    mass_[s] += masses[k];
  }

  if (dc_dx != nullptr)
  {
    // compute the and center of mass of each voronoi cell
    centroid_.clear();
    centroid_.resize( delaunay_.nb()*dim , 0.0 );
    for (index_t k = 0; k < simplices_.nb(); k++)
    {
      index_t s  = simplices_.simplex2site(k);
      for (coord_t d = 0; d < dim; d++)
        centroid_[ s*dim + d ] += moments[ k*dim + d];
    }

    // compute the centroid and gradient
    for (index_t k = 0; k < delaunay_.nb(); k++)
    for (coord_t d = 0; d < dim; d++)
    {
      // if weights are too high, some cells could be submerged which gives a zero volume
      if (mass_[k] < 1e-12) centroid_[k*dim+d] = delaunay_[k][d];
      else
        centroid_[k*dim+d] /= mass_[k];
      dc_dx[ k*dim + d] = 2.0*mass_[k]*( delaunay_[k][d] - centroid_[k*dim+d] );
    }
  }

  if (dc_dw != nullptr)
  {
    // compute the gradient of the transport map
    for (index_t k = 0; k < delaunay_.nb(); k++)
      dc_dw[k] = nu_[k] - mass_[k];
  }

  // get the energy
  t0 = clock();
  typedef voronoi::Integrand_Transport_Energy Integrand_t;
  Integrand_t integrand(delaunay_,simplices_,*density_,dim);
  Functional<Integrand_t> f(integrand,simplices_);
  f.integrate(simplices_);
  real_t energy = f.value();
  time_integrate += real_t( clock() - t0 )/real_t(CLOCKS_PER_SEC);

  // add te contribution from the weights and dirac masses
  for (index_t k=0;k<delaunay_.nb();k++)
    energy += weight_[k]*(nu_[k] - mass_[k]);

  // calculate the norm of the gradient and report convergence criteria
  real_t gnorm;
  if (mode_ == 0)
    gnorm = calculate_norm(dc_dx,delaunay_.nb()*dim);
  else
    gnorm = calculate_norm(dc_dw,delaunay_.nb());


  time_integrate /= ProcessCPU::maximum_concurrent_threads();
  real_t mass_min = * std::min_element( mass_.begin() , mass_.end() );
  real_t mass_max = * std::max_element( mass_.begin() , mass_.end() );

  if (print_)
  {
    printf("%6lu|%10.3e|%10.3e|%11.1e|%11.1e|%11.1e|%11.1e|%11.3e|%11.3e|%6lu|%6lu|%11.1e\n",
      iteration_,
      energy,
      gnorm,
      diagram_.time_voronoi(),
      diagram_.time_neighbours(),
      diagram_.time_decompose(),
      time_integrate,
      mass_min,
      mass_max,
      diagram_.minimum_neighbours(),
      diagram_.maximum_neighbours(),
      diagram_.average_neighbours()
    );
  }

  properties_.mass_min.push_back( mass_min );
  properties_.mass_max.push_back( mass_max );
  properties_.time_voronoi.push_back( diagram_.time_voronoi() );
  properties_.time_neighbours.push_back( diagram_.time_neighbours() );
  properties_.time_integration.push_back( time_integrate );
  properties_.time_triangulation.push_back( diagram_.time_decompose() );
  properties_.neighbours_average.push_back( diagram_.average_neighbours() );
  properties_.neighbours_minimum.push_back( diagram_.minimum_neighbours() );
  properties_.neighbours_maximum.push_back( diagram_.maximum_neighbours() );
  properties_.energy.push_back( energy );
  properties_.gradient.push_back( gnorm );
  time_integrate_ = time_integrate;

  if (iteration_ % save_every_ == 0) save_snapshot();

  return energy;
}

void
SDOT_Properties::save( const std::string& filename ) const
{
  json J;
  J["mass_min"] = mass_min;
  J["mass_max"] = mass_max;
  J["time_voronoi"] = time_voronoi;
  J["time_neighbours"] = time_neighbours;
  J["integration"] = time_integration;
  J["triangulation"] = time_integration;
  J["neighbours_average"] = neighbours_average;
  J["neighbours_minimum"] = neighbours_minimum;
  J["neighbours_maximum"] = neighbours_maximum;
  J["energy"] = energy;
  J["gradient"] = gradient;
  std::ofstream output(filename);
  output << std::setw(4) << J << std::endl;
  output.close();
}

void
SDOT_Snapshot::save( const std::string& filename ) const
{
  json J;
  J["dim"] = dim;
  J["points"] = points;
  J["weights"] = weights;
  J["neighbours"] = neighbours;
  J["density"] = name;
  J["mass"] = mass;
  std::ofstream output(filename);
  output << std::setw(4) << J << std::endl;
  output.close();
}

template<typename type>
void
SemiDiscreteOptimalTransport<type>::save_snapshot()
{
  snapshot_.clear();
  std::vector<real_t> coordinates(delaunay_.nb() * domain_.number());
  for (index_t k = 0, j = 0 ; k < delaunay_.nb(); k++)
  for (index_t d = 0; d < domain_.number(); d++)
    coordinates[j++] = delaunay_[k][d];
  snapshot_.dim = domain_.number();
  snapshot_.points = coordinates;
  snapshot_.weights = weight_;
  snapshot_.neighbours = diagram_.neighbour_counts();
  snapshot_.mass = mass_;
  snapshot_.name = diagram_.name();
  snapshot_.save(prefix_ + "-snapshot-iter-" + std::to_string(iteration_) + ".json" );
}

template<typename type>
struct nlopt_data_optimal_transport
{
	SemiDiscreteOptimalTransport<type>& transport;
	index_t eval_count;
	real_t objective;
  int mode; // 0 for coordinates , 1 for weights
};

template<typename type>
void
SemiDiscreteOptimalTransport<type>::set_delaunay( const real_t* x , coord_t dim )
{
  for (index_t k=0;k<delaunay_.nb();k++)
  for (coord_t d=0;d<dim;d++)
    delaunay_[k][d] = x[k*dim+d];
}

template<typename type>
void
SemiDiscreteOptimalTransport<type>::set_weights( const real_t* w )
{
  for (index_t k=0;k<delaunay_.nb();k++)
    weight_[k] = w[k];
  coord_t dim = delaunay_.dim() -1;
  real_t wm = * std::max_element( weight_.begin() , weight_.end() );
  //printf("setting dim = %u, wm = %g\n",dim,wm);
  for (index_t k=0;k<weight_.size();k++)
    delaunay_[k][dim] = std::sqrt( std::max( wm - weight_[k] , 0.0 ) );
}

template<typename type>
real_t
SemiDiscreteOptimalTransport<type>::transport_objective( index_t n , const double* x , double* grad )
{
  // set the coordinates of the delaunay points
  if (mode_ == 0)
    set_delaunay(x,domain_.number());
  else if (mode_ == 1)
    set_weights(x);
  else
    avro_assert_not_reached;

  // recompute the voronoi diagram
  compute_laguerre();

  index_t nb_points = delaunay_.nb();
  coord_t dim = domain_.number();

  // compute the energy and gradients
  std::vector<real_t> dE_dZ(nb_points*dim),dE_dW(nb_points);
  real_t energy = evaluate( dE_dZ.data() , dE_dW.data() );
  if (grad != nullptr)
  {
    if (mode_ == 0)
    {
      // cvt mode
      for (index_t k=0;k<dE_dZ.size();k++)
        grad[k] = dE_dZ[k];
    }
    else if (mode_ == 1)
    {
      // otm mode: flip the sign to make it convex
      energy *= -1.0;
      for (index_t k=0;k<dE_dW.size();k++)
        grad[k] = -dE_dW[k];
    }
    else
      avro_assert_not_reached;
  }
  return energy;
}

template<typename type>
double
nlopt_transport_objective( unsigned n , const double* x , double* grad, void* data0 )
{
  nlopt_data_optimal_transport<type>* data = (nlopt_data_optimal_transport<type>*)(data0);
  SemiDiscreteOptimalTransport<type>& transport = data->transport;
  transport.iteration()++;
  return transport.transport_objective( n , x , grad );
}

#define OPTIMIZER_NLOPT 1
OptimalTransportBase* __transport__ = nullptr;

template<typename type>
void
SemiDiscreteOptimalTransport<type>::start()
{
  printf("--------------------------------------------------------------------------------------\n");
  printf("%5s|%10s|%10s|%10s|%10s|%10s|%10s|%10s\n","iter  "," energy "," gradient "," t_vor (s) "," t_knn (s) " , " t_tri (s) " , " t_int (s) ", " min(m) ");
  printf("--------------------------------------------------------------------------------------\n");
  iteration_ = 0;
}

template<typename type>
void
SemiDiscreteOptimalTransport<type>::optimize_points( index_t nb_iter )
{
  properties_.clear();

  mode_ = 0;
  coord_t dim = domain_.number();
  index_t n = delaunay_.nb()*dim;

  std::vector<real_t> x(n);
  for (index_t k=0;k<delaunay_.nb();k++)
  for (index_t d=0;d<dim;d++)
    x[k*dim+d] = delaunay_(k,d);

  // setup the optimizer
  nlopt::opt opt( nlopt::LD_LBFGS , n );
  nlopt_data_optimal_transport<type> data = {*this,0,1,0};

	// set the objective function
	opt.set_min_objective( &nlopt_transport_objective<type> , static_cast<void*>(&data) );

  // set some optimization parameters
  opt.set_xtol_rel(1e-12);
  opt.set_ftol_rel(1e-12);
  opt.set_maxeval(nb_iter);

  // set the lower and upper bounds on the weights
  std::vector<real_t> lower_bound( n , 0.0 );
  std::vector<real_t> upper_bound( n ,  1.0 );
  opt.set_lower_bounds(lower_bound);
  opt.set_upper_bounds(upper_bound);

  real_t f_opt;
  nlopt::result result = nlopt::result::SUCCESS;
  start();
  result = opt.optimize(x,f_opt);
  std::string desc = nloptResultDescription(result);
  printf("nlopt result: %s\n",desc.c_str());

  save_snapshot();
}

template<typename type>
void
SemiDiscreteOptimalTransport<type>::optimize_points_lloyd( index_t nb_iter )
{
  properties_.clear();

  mode_ = 0;
  coord_t dim = domain_.number();
  index_t n = delaunay_.nb()*dim;

  centroid_.clear();
  centroid_.resize(n);
  for (index_t k=0;k<delaunay_.nb();k++)
  for (index_t d=0;d<dim;d++)
    centroid_[k*dim+d] = delaunay_(k,d);

  start();
  for (iteration_ = 0; iteration_ < nb_iter ; iteration_++)
  {
    // set the current delaunay vertices
    set_delaunay( centroid_.data() , dim );

    // evaluate the voronoi diagram for these points
    std::vector<real_t> dc_dx(n); // not used
    evaluate( dc_dx.data() , nullptr );
  }

  save_snapshot();
}

template<typename type>
void
SemiDiscreteOptimalTransport<type>::optimize_weights( index_t nb_iter )
{
  properties_.clear();

  mode_ = 1;
  index_t n = delaunay_.nb();
  std::vector<real_t> x(n);
  for (index_t k=0;k<delaunay_.nb();k++)
    x[k] = 0.0;//mass_[k];//weight_[k];

  // setup the optimizer
  nlopt::opt opt( nlopt::LD_LBFGS , n );
  nlopt_data_optimal_transport<type> data = {*this,0,1,1};

	// set the objective function
	opt.set_min_objective( &nlopt_transport_objective<type> , static_cast<void*>(&data) );

  // set some optimization parameters
  opt.set_xtol_rel(1e-16);
  opt.set_ftol_rel(1e-16);
  opt.set_maxeval(nb_iter);

  // set the lower and upper bounds on the weights
  std::vector<real_t> lower_bound( n , 0.0 );
  std::vector<real_t> upper_bound( n ,  weight_max_ );
  opt.set_lower_bounds(lower_bound);
  opt.set_upper_bounds(upper_bound);

  real_t f_opt;
  nlopt::result result = nlopt::result::SUCCESS;
  start();
  result = opt.optimize(x,f_opt);
  std::string desc = nloptResultDescription(result);
  printf("nlopt result: %s\n",desc.c_str());

  save_snapshot();
}

template<typename type>
void
SemiDiscreteOptimalTransport<type>::generate_bluenoise()
{
  // optimize the points
  optimize_points(100);

  // compute the target mass
  real_t m_total = 0.0;
  for (index_t k = 0; k < mass_.size(); k++)
    m_total += mass_[k];
  nu_.clear();
  nu_.resize( delaunay_.nb() , m_total/delaunay_.nb() );

  printf("total mass = %g, target mass = %g\n",m_total,m_total/delaunay_.nb());

  // generate blue noise
  for (index_t iter = 0; iter < 5; iter++)
  {
    // compute the transport map
    //optimize_weights(10);

    // optimize the points
    //optimize_points_lloyd(1);
  }

  // generate blue noise
  for (index_t iter = 0; iter < 5; iter++)
  {
    // compute the transport map
    //optimize_weights(10);

    // optimize the points
  //  optimize_points(2);
  }
}

template<typename type>
void
SemiDiscreteOptimalTransport<type>::stochastic_gradient_descent( index_t nb_iter )
{

  start();
  mode_ = 0;

  index_t n = delaunay_.nb();
  coord_t dim = domain_.number();

  // initialize the gradient
  std::vector<real_t> f(n, 0.0);
  std::vector<real_t> dE_dw( n ,0.0 );
  std::vector<real_t> a( n , 1.0 );
  for (iteration_ = 0; iteration_ < nb_iter; iteration_++)
  {
    // sample the density measure
    std::vector<real_t> y(dim,0.0);

    // determine which laguerre cell contains y
    index_t i = 0;

    // update the gradient
    dE_dw = a;
    dE_dw[i] -= 1.0;

    real_t l0 = 10.0;
    real_t tau = 0.1/( 1.0 + iteration_/l0 );

    f[i] += tau*dE_dw[i];

    // set the weights and re-evaluate
    set_weights( a.data() );
    evaluate( nullptr , nullptr );
  }

}

template<typename type>
std::vector<real_t>
SemiDiscreteOptimalTransport<type>::get_sites() const {
  index_t nb_points = delaunay_.nb();
  index_t dim = domain_.number(); // maybe make dim a parameter to this function
  std::vector<real_t> sites( nb_points*dim );
  index_t i = 0;
  for (index_t k = 0; k < nb_points; k++)
  for (coord_t d = 0; d < dim; d++)
    sites[i++] = delaunay_[k][d];
  return sites;
}

template class LaguerreDiagram<Simplex>;
template class LaguerreDiagram<Polytope>;
template class SemiDiscreteOptimalTransport<Simplex>;
template class SemiDiscreteOptimalTransport<Polytope>;

} // voronoi

} // avro
