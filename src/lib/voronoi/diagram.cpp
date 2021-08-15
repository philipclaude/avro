#include "numerics/geometry.h"
#include "numerics/nlopt_result.h"

#include "voronoi/diagram.h"

#include <nlopt.hpp>

#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/IterativeLinearSolvers>

#include <geogram/nl.h>

namespace avro
{

namespace voronoi
{

#define PARALLEL 1

void
PowerDiagram::initialize() {
  // creates voronoi cells and sites the target mass and weight arrays
  // should only be called once
  cell_.resize( sites_.nb() );
  for (index_t k = 0; k < sites_.nb(); k++) {
    cell_[k] = std::make_shared<Cell>( k , sites_ , domain_ , *search_ );
    cell_[k]->set_ambient_dimension(ambient_dim_);
  }
  nu_.resize(sites_.nb() , 0.0 );
  weight_.resize( sites_.nb() , 0.0 );
}

void
PowerDiagram::allocate_sites( index_t n ) {
  sites_.clear();
  std::vector<real_t> x(sites_.dim());
  for (index_t k = 0; k < n; k++)
    sites_.create(x.data());
}

void
PowerDiagram::set_sites( const real_t* x ) {

  // copy the points
  index_t i = 0;
  for (index_t k = 0; k < sites_.nb(); k++)
  for (coord_t d = 0; d < ambient_dim_; d++)
    sites_[k][d] = x[i++];

  // set the points into the nearest neighbour search structure
  if (search_ == nullptr)
    search_ = GEO::NearestNeighborSearch::create(sites_.dim(),"ANN");
  search_->set_points( sites_.nb() , sites_[0] );
}

void
PowerDiagram::set_sites( const Points& p ) {

  coord_t dim = p.dim();
  avro_assert( dim == sites_.dim() );

  // copy the points and set the coordinates into the search structure
  p.copy(sites_);
  if (search_ == nullptr)
    search_ = GEO::NearestNeighborSearch::create(dim,"ANN");
  search_->set_points( sites_.nb() , sites_[0] );
}

void
PowerDiagram::compute() {

  // clear any previous power diagram data
  Topology<Polytope>::clear();
  vertices_.clear();
  triangles_.clear();
  triangle2site_.clear();
  edges_.clear();
  polytope2site_.clear();

  // clear any optimization data
  de_dx_.resize( sites_.nb() * ambient_dim_ , 0.0 );
  de_dw_.resize( sites_.nb() , 0.0 );
  centroid_.resize( sites_.nb() * ambient_dim_ , 0.0 );
  cell_volume_.resize( sites_.nb() , 0.0 );

  // reset the energy, volume and area
  energy_ = 0.0;
  volume_ = 0.0;
  area_   = 0.0;

  // compute the power diagram
  //printf("computing the power diagram in %ud\n",sites_.dim());
  clock_t t0 = clock();
  #if PARALLEL
  typedef PowerDiagram thisclass;
  ProcessCPU::parallel_for(
    parallel_for_member_callback( this , &thisclass::compute ), 0,cell_.size()
  );
  #else
  for (index_t k = 0; k < cell_.size(); k++) {
    compute(k);
  }
  #endif

  // record the time it took to compute the voronoi diagram
  time_voronoi_ = real_t(clock()-t0)/real_t(CLOCKS_PER_SEC)/ProcessCPU::maximum_concurrent_threads();
}


void
PowerDiagram::compute( index_t k ) {

  // get the candidate list of elements in the background mesh to clip against
  // (only one will be used)
  std::vector<index_t> ball;
  domain_.get_candidate_elems( sites_[k] , ball );

  // clip the voronoi cell against the mesh
  cell_[k]->compute( ball );

  // retrieve the mass and moment data
  const std::vector<real_t>& moment = cell_[k]->moment();
  real_t volume_k = cell_[k]->volume();
  cell_volume_[k] = volume_k;

  #if 0 // this check should only be used when weights are 0, oterwise cells could indeed vanish
  if (volume_k <= 0.0) {
    sites_.print(k);
  }
  avro_assert( volume_k > 0.0 );
  #endif

  // accumulate the total volume and boundary area (used for testing)
  volume_ += volume_k;
  area_   += cell_[k]->boundary_area();

  // compute the centroid of the cell
  if (volume_k > 0.0) {
    for (coord_t d = 0; d < ambient_dim_; d++)
      centroid_[k*ambient_dim_+d] = moment[d]/volume_k;
  }
  else {
    for (coord_t d = 0; d < ambient_dim_; d++)
      centroid_[k*ambient_dim_+d] = sites_[k][d]; // retain the previous value
  }

  // compute the gradient of the energy with respect to the site locations
  for (coord_t d = 0; d < ambient_dim_; d++)
    de_dx_[k*ambient_dim_+d] = 2.0*( volume_k*sites_[k][d] - moment[d] );

  // compute the gradient of the energy with respect to the weights
  de_dw_[k] = nu_[k] - volume_k;

  // add the contribution of the cell energy to the total
  energy_   += cell_[k]->energy() + weight_[k] * (nu_[k] - volume_k );
}

void
PowerDiagram::accumulate(bool complete) {
  // accumulate all the cells into a polytope mesh that we can visualize
  for (index_t k = 0; k < cell_.size(); k++) {
    add_cell( *cell_[k].get() , complete );
    cell_[k]->clear();
  }
}

void
PowerDiagram::add_cell( const voronoi::Cell& cell , bool complete ) {

  // add the points
  index_t nb_points = vertices_.nb();
  for (index_t j = 0; j < cell.points().nb(); j++) {
    vertices_.create( cell.points()[j] );
    if (complete) {
      std::vector<int> b = cell.points().incidence().get(j);
      vertices_.incidence().add( b.data() , b.size() );
    }
  }

  // add a dummy polytope
  // (we don't actually need the polytope since we already have the triangles to visualize)
  for (index_t j = 0; j < cell.nb(); j++) {
    if (complete) {
      std::vector<index_t> polytope = cell.get(j);
      for (index_t i = 0; i < polytope.size(); i++)
        polytope[j] += nb_points;
        add( polytope.data() , polytope.size() );
    }
    else {
      index_t dummy = 0;
      add( &dummy , 1 );
    }
    polytope2site_.push_back(cell.site());
  }

  // add the triangle data
  const std::vector<index_t>& t = cell.triangles();
  for (index_t j = 0; j < t.size(); j++)
    triangles_.push_back(t[j]+nb_points);
  for (index_t j = 0; j < t.size()/3; j++)
    triangle2site_.push_back(cell.site());

  // add the edge data
  const std::vector<index_t>& e = cell.edges();
  for (index_t j = 0; j < e.size(); j++)
    edges_.push_back(e[j]+nb_points);
}

void
PowerDiagram::create_field() {
  // creates a field so that we can visualize each voronoi cell with a constant colour
  site_field_ = std::make_shared<SiteField>(*this);
  fields().make("sites",site_field_);
}

struct nlopt_power_diagram_data {
	PowerDiagram& diagram;
	index_t eval_count;
	real_t objective;
  int mode; // 0 for coordinates , 1 for weights
};

double
nlopt_objective( unsigned n , const double* x , double* grad, void* data0 )
{
  nlopt_power_diagram_data* data = (nlopt_power_diagram_data*)(data0);
  PowerDiagram& diagram = data->diagram;
  diagram.iteration()++;
  return diagram.evaluate_objective( n , x , grad );
}

void
PowerDiagram::start() {
  if (verbose_) {
    printf("-------------------------------------------------\n");
    printf("%8s|%12s|%12s|%12s\n","iter  "," energy "," gradient "," time (s)");
    printf("-------------------------------------------------\n");
  }
  iteration_ = 0;
  sub_iteration_ = 0;
}

static real_t
calculate_norm( const real_t* x , index_t nb ) {
  real_t n = 0.0;
  for (index_t k = 0; k < nb; k++)
    n += x[k]*x[k];
  return std::sqrt(n);
}

real_t
PowerDiagram::evaluate_objective( index_t n , const real_t* x , real_t* grad ) {

  real_t gnorm = 0.0;
  if (mode_ == 0) {
    set_sites(x);
  }
	if (mode_ == 1) {
		for (index_t k = 0; k < n; k++)
			weight_[k] = x[k];

		// lift the points
		coord_t lifted_dim = sites_.dim()-1;
		real_t wmax = * std::max_element( weight_.begin() , weight_.end() );
		for (index_t k = 0; k < sites_.nb(); k++)
			sites_[k][lifted_dim] = std::sqrt( wmax - weight_[k] );
		if (search_ == nullptr)
			search_ = GEO::NearestNeighborSearch::create(sites_.dim(),"ANN");
		search_->set_points( sites_.nb() , sites_[0] );
	}

  compute();

	real_t regularization_total = 0.0;
	real_t multiplier = 1.0;
	if (grad != nullptr) {
	  if (mode_ == 0) {
	    for (index_t i = 0; i < de_dx_.size(); i++)
	      grad[i] = de_dx_[i];
	  }
		else if (mode_ == 1) {
			multiplier = -1.0;
			for (index_t i = 0; i < de_dw_.size(); i++)
				grad[i] = multiplier*de_dw_[i];
		}
	}

	if (mode_ == 1) {
		// add a regularization term so that the hessian is not singular
		real_t epsilon = regularization();
		for (index_t k = 0; k < n; k++) {
			real_t dg = epsilon * nu_[k] * weight_[k];
			regularization_total += dg;
			de_dw_[k] += dg;
		}
	}

  gnorm = 1.0;
	if      (mode_ == 0) gnorm = calculate_norm(de_dx_.data() , de_dx_.size() );
	else if (mode_ == 1) gnorm = calculate_norm(de_dw_.data() , de_dw_.size() );

	real_t min_vol = *std::min_element( cell_volume_.begin() , cell_volume_.end() );
	real_t max_vol = *std::max_element( cell_volume_.begin() , cell_volume_.end() );

  if (verbose_)
  printf("%5lu-%3lu|%12.3e|%12.3e|%12.3e|%12.3e|%12.3e|%12.3e\n",
    iteration_,sub_iteration_,
    energy_,
    gnorm,
    time_voronoi_,
		min_vol,max_vol,regularization_total
  );

  return multiplier*energy_;
}

void
PowerDiagram::optimize_points( index_t nb_iter ) {

  mode_ = 0;
  coord_t dim = ambient_dim_;
  index_t n = sites_.nb()*dim;

  std::vector<real_t> x(n);
  for (index_t k = 0; k < sites_.nb(); k++)
  for (index_t d = 0; d < dim; d++)
    x[k*dim+d] = sites_(k,d);

  // setup the optimizer
  nlopt::opt opt( nlopt::LD_LBFGS , n );
  nlopt_power_diagram_data data = {*this,0,1,0};

  // set the objective function
  opt.set_min_objective( &nlopt_objective , static_cast<void*>(&data) );

  // set some optimization parameters
  opt.set_xtol_rel(1e-12);
  opt.set_ftol_rel(1e-12);
  opt.set_maxeval(nb_iter);

  // set the lower and upper bounds on the coordinates
  #if 0 // this should not be used for general meshes
  std::vector<real_t> lower_bound( n , 0.0 );
  std::vector<real_t> upper_bound( n ,  1.0 );
  opt.set_lower_bounds(lower_bound);
  opt.set_upper_bounds(upper_bound);
  #endif

  real_t f_opt;
  nlopt::result result = nlopt::result::SUCCESS;
  start();
  try {
    result = opt.optimize(x,f_opt);
    std::string desc = nloptResultDescription(result);
    printf("nlopt result: %s\n",desc.c_str());
  }
	catch (std::exception& e) {
		std::cout << e.what() << std::endl;
		printf("runtime error\n");
	}

}

void
PowerDiagram::optimize_points_lloyd( index_t nb_iter ) {

	mode_ = 0;
  coord_t dim = ambient_dim_;
  index_t n = sites_.nb()*dim;

	// set the initial sites into the design variables
  std::vector<real_t> x(n);
  for (index_t k = 0; k < sites_.nb(); k++)
  for (index_t d = 0; d < dim; d++)
    x[k*dim+d] = sites_(k,d);

	start();
	for (index_t k = 0; k < nb_iter; k++) {

			evaluate_objective( n , x.data() , nullptr );
			for (index_t i = 0; i < n; i++)
				x[i] = centroid_[i];
			iteration_++;
	}
}

void
PowerDiagram::optimize_weights( index_t nb_iter , const std::vector<real_t>& mass ) {

	mode_ = 1;
  index_t n = sites_.nb();

	// set the initial weights to zero and save the target mass
  std::vector<real_t> x(n,0.0);
	for (index_t k = 0; k < n; k++)
		nu_[k] = mass[k];

  // setup the optimizer
  nlopt::opt opt( nlopt::LD_LBFGS , n );
  nlopt_power_diagram_data data = {*this,0,1,int(mode_)};

  // set the objective function
  opt.set_min_objective( &nlopt_objective , static_cast<void*>(&data) );

  // set some optimization parameters
  opt.set_xtol_rel(1e-12);
  opt.set_ftol_rel(1e-12);
  opt.set_maxeval(nb_iter);

  // set the lower and upper bounds on the weights
  #if 1
  std::vector<real_t> lower_bound( n , 0.0 );
	opt.set_lower_bounds(lower_bound);
  //std::vector<real_t> upper_bound( n ,  1e6 );
  //opt.set_upper_bounds(upper_bound);
  #endif

  real_t f_opt;
  nlopt::result result = nlopt::result::SUCCESS;
  start();
  try {
    result = opt.optimize(x,f_opt);
    std::string desc = nloptResultDescription(result);
    printf("nlopt result: %s\n",desc.c_str());
  }
	catch (std::exception& e) {
		std::cout << e.what() << std::endl;
		printf("runtime error\n");
	}
}

std::string
eigen_error_message( const Eigen::ComputationInfo& info ) {
  if (info == Eigen::Success) return "success";
  if (info == Eigen::NumericalIssue) return "numerical issue";
  if (info == Eigen::InvalidInput) return "invalid input";
  if (info == Eigen::NoConvergence) return "no convergence";
  return "unknown error";
}

class LinearSolver_Eigen {

public:

	void initialize( index_t m , index_t n , index_t nnz ) {
		matrix_.setZero();
		matrix_.resize(m,n);
		matrix_.reserve(nnz);
	}

	void fill( const std::vector< std::tuple<index_t,index_t,real_t>>& values ) {

		typedef Eigen::Triplet<real_t> T;
		std::vector<T> triplets(values.size());
		for (index_t k = 0; k < values.size(); k++) {
			index_t i = std::get<0>(values[k]);
			index_t j = std::get<1>(values[k]);
			real_t  v = std::get<2>(values[k]);
			triplets[k] = T(i,j,v);
		}
		matrix_.setFromTriplets(triplets.begin(),triplets.end());
	}

	void begin_matrix() {}
	void end_matrix() {}

	void set_rhs( const std::vector<real_t>& b ) {
		rhs_ = b;
	}

	void add_entry( index_t i , index_t j , real_t value ) {
		matrix_.coeffRef(i,j) = value;
	}

	void solve( std::vector<real_t>& x ) {

		// set the right-hand-side vector
		index_t n = rhs_.size();
		Eigen::VectorXd b(n);
		for (index_t k = 0; k < n; k++)
			b(k) = rhs_[k];

    //matrix_.makeCompressed();
		Eigen::SparseLU<Eigen::SparseMatrix<real_t> > solver;
		//Eigen::BiCGSTAB<Eigen::SparseMatrix<real_t> , Eigen::DiagonalPreconditioner<real_t> > solver; //Eigen::IncompleteLUT<real_t> > solver;
		//Eigen::ConjugateGradient<Eigen::SparseMatrix<real_t> , Eigen::Lower | Eigen::Upper , Eigen::DiagonalPreconditioner<real_t> > solver;
    //Eigen::SparseQR<Eigen::SparseMatrix<real_t> , Eigen::COLAMDOrdering<int> > solver;
    //Eigen::LeastSquaresConjugateGradient< Eigen::SparseMatrix<real_t> , Eigen::DiagonalPreconditioner<real_t> > solver;

		// set some parameters
    //solver.setTolerance(1e-3);
    //solver.setMaxIterations(10000);

		// factorize/precondition the matrix
		//clock_t t0 = clock();
		solver.compute(matrix_);
		if (solver.info() != Eigen::Success) {
			printf("decomposition failed\n");
      avro_assert_not_reached;
		}

		// solve the system
		Eigen::VectorXd dw(n);
		dw = solver.solve(b);
		if (solver.info() != Eigen::Success) {
			printf("solve failed\n");
      std::cout << eigen_error_message(solver.info()) << std::endl;
      //printf("iterations = %lu, error = %g",solver.iterations(),solver.error());
      avro_assert_not_reached;
		}
		//clock_t t1 = clock();
		//printf("--> linear solve time = %g sec.\n",real_t(t1-t0)/real_t(CLOCKS_PER_SEC));

		// save the result
		for (index_t k = 0; k < n; k++)
			x[k] = dw(k);
	}

private:
	Eigen::SparseMatrix<real_t> matrix_;
	std::vector<real_t> rhs_;
};

#if 1
class LinearSolver_OpenNL {

public:

	void initialize( index_t m , index_t n , index_t nnz ) {
		avro_assert(m == n); // for now
		nlNewContext();

		nlSolverParameteri( NL_NB_VARIABLES,NLint(n) );
		nlSolverParameteri( NL_SOLVER , NL_CG );
		nlSolverParameteri( NL_PRECONDITIONER , NL_PRECOND_JACOBI );
		nlSolverParameteri( NL_SYMMETRIC , NL_TRUE );
		//nlSolverParameteri( NL_THRESHOLD , 1e-5 );
		//nlSolverParameteri( NL_MAX_ITERATIONS , 10000 );

		//nlEnable( NL_VARIABLES_BUFFER );
		//nlEnable( NL_NO_VARIABLES_INDIRECTION );

		nlBegin(NL_SYSTEM);
	}

	void begin_matrix() {
		nlBegin(NL_MATRIX);
	}

	void end_matrix() {
		nlEnd(NL_MATRIX);
	}

	void set_rhs( const std::vector<real_t>& b ) {

		//nlBindBuffer(NL_VARIABLES_BUFFER,0,x_.data(),NLuint(sizeof(real_t)));
		for (index_t k = 0; k < b.size(); k++)
			nlAddIRightHandSide(k,b[k]);
		printf("added rhs\n");
	}

	void add_entry( index_t i , index_t j , real_t value ) {
		nlAddIJCoefficient(i,j,value);
	}

	void solve( std::vector<real_t>& x ) {

		nlEnd(NL_SYSTEM);
		nlSolve();

		index_t n = x.size();
		for (index_t k = 0; k < n; k++)
			x[k] = nlGetVariable(k);

		print_inline(x);

		nlDeleteContext(nlGetCurrent());

	}

private:
	std::vector<real_t> x_;

};
#endif

template<typename LinearSolver>
class NewtonStep {

public:
	NewtonStep( const PowerDiagram& diagram ) :
		diagram_(diagram),
		regularization_(diagram.regularization())
	{}

	void set_gradient( const std::vector<real_t>& grad ) {
		gradient_ = grad;
	}

	void set_regularization( real_t x ) {
		regularization_ = x;
	}

	void get_pattern( std::vector< std::vector<index_t> >& pattern ) {

		index_t n = diagram_.sites().nb();

		// get the list of all i,j entries
		std::vector< std::vector<index_t> > rows( n );

		for (index_t k = 0; k < n; k++) {

			rows[k].reserve(20*k);
			rows[k].push_back(k);

			// get the power facets
			const std::map<int,PowerFacet>& facets = diagram_.cell(k).facets();
			for (std::map<int,PowerFacet>::const_iterator it = facets.begin(); it != facets.end(); ++it) {

				const PowerFacet& facet = it->second;
				if (facet.pj < 0) continue; // skip boundary facets

				int pi = facet.pi;
				int pj = facet.pj;
				avro_assert_msg( pi == int(k) , "pi = %d, k = %lu" , pi , k );
				avro_assert_msg( pi != pj , "pi = %d, pj = %d" , pi , pj );

				rows[k].push_back(pj);
			}
		}
	}

	void build_hessian() {

		index_t n = diagram_.sites().nb();
		coord_t dim = diagram_.ambient_dimension();

		// allocate the hessian matrix
		linear_solver_.initialize( n , n , 20*n );

		std::vector<std::vector<index_t>> pattern;
		get_pattern(pattern);

		linear_solver_.begin_matrix();
		linear_solver_.set_rhs( gradient_ );

		// loop through the power cells
		std::vector< std::tuple<index_t,index_t,real_t> > triplets;
		triplets.reserve(20*n);

		for (index_t k = 0; k < n; k++) {

			const real_t* xi = diagram_.sites()[k];
			real_t Hii = 0.0;

			// get the voronoi cell
			const Cell& cell = diagram_.cell(k);

			// get the power facets
			const std::map<int,PowerFacet>& facets = cell.facets();

			for (std::map<int,PowerFacet>::const_iterator it = facets.begin(); it != facets.end(); ++it) {

				const PowerFacet& facet = it->second;
				if (facet.pj < 0) continue; // skip boundary facets

				int pi = facet.pi;
				int pj = facet.pj;
				avro_assert_msg( pi == int(k) , "pi = %d, k = %lu" , pi , k );
				avro_assert_msg( pi != pj , "pi = %d, pj = %d" , pi , pj );

				const real_t* xj = diagram_.sites()[pj];
				real_t Hij = 0.5*facet.Aij / numerics::distance(xi,xj,dim);

				// set the off-diagonal term
				triplets.push_back( std::tuple<index_t,index_t,real_t>(k,pj,Hij) );

				// add the contribution to the diagonal term
				Hii -= Hij;
			}

			// add a regularization term
			real_t dh = regularization_ * diagram_.nu(k);
			Hii += dh;

			// set the diagonal term
			triplets.push_back( std::tuple<index_t,index_t,real_t>(k,k,Hii) );
		}
		linear_solver_.end_matrix();
		linear_solver_.fill(triplets);
	}

	void solve( std::vector<real_t>& pk ) {

		index_t n = diagram_.sites().nb();
		linear_solver_.solve(pk);
		for (index_t k = 0; k < n; k++)
			pk[k] *= -1.0;
	}

private:
	const PowerDiagram& diagram_;
	LinearSolver linear_solver_;
	std::vector<real_t> gradient_;

	real_t regularization_;
};

bool
PowerDiagram::optimize_weights_kmt( index_t nb_iter , const std::vector<real_t>& mass ) {

	mode_ = 1; // we are in site-optimization mode
  index_t n = sites_.nb(); // number of design variables

	// set the initial weights to zero and save the target mass
  std::vector<real_t> x(n,0.0);
	for (index_t k = 0; k < n; k++)
		nu_[k] = mass[k];

  // print the optimization header
	start();

  // compute the voronoi diagram (weights = 0) which also evaluates the gradient
  evaluate_objective( n , x.data() , nullptr );

  // determine the a0 constant in the Kitagawa algorithm
  real_t vmin = *std::min_element( cell_volume_.begin() , cell_volume_.end() );
  real_t mmin = *std::min_element( nu_.begin() , nu_.end() );
  real_t gnorm = -1.0, gnorm0 = -1.0;
  real_t a0 = 0.5*( std::min(vmin,mmin) );
  if (verbose_) printf("--> determined a0 = %g\n",a0);

  bool converged = true;
	for (index_t iter = 0; iter < nb_iter; iter++) {

		// build the hessian matrix
		NewtonStep<LinearSolver_Eigen> newton(*this);
		//NewtonStep<LinearSolver_OpenNL> newton(*this);
		newton.set_gradient(de_dw_);
		newton.build_hessian();

		// determine the descent direction
		std::vector<real_t> pk( n );
		newton.solve(pk);

    // determine the step size
    real_t alpha = 1.0;

    // copy the original weights, and calculate the current norm of the energy gradient
    std::vector<real_t> x0(x.begin(),x.end());
    gnorm0 = calculate_norm( de_dw_.data() , n );
		if (gnorm0 < 1e-16) {
      gnorm = gnorm0;
      break;
    }

    sub_iteration_ = 0;
    while (true) {

      sub_iteration_++;

      // update the weights
  		for (index_t k = 0; k < n; k++)
  			x[k] = x0[k] + alpha*pk[k];

      // re-evaluate at these weights
  		evaluate_objective( n , x.data() , nullptr );

      // check if we have sufficient decrease in the gradient
      // as well as the minimum volume isn't too small
      gnorm = calculate_norm( de_dw_.data() , n );
      vmin  = *std::min_element( cell_volume_.begin() , cell_volume_.end() );

      // check for sufficient decrease in the gradient, and the volume isn't too small
      if (gnorm <= (1.0 - alpha/2.0)*gnorm0 && vmin > a0) break;

      // try halving the step size
      alpha = alpha / 2.0;

      if (alpha < 1e-7 || sub_iteration_ > 10) {
				converged = false;
				break;
			}
    }
		if (!converged) break;
    if (gnorm < 1e-10) {
      converged = true;
      break;
    }

		iteration_++;
	}

  residual_ = gnorm;

  return converged;

}

} // voronoi

} // avro
