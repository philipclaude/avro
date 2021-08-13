#include "numerics/geometry.h"
#include "numerics/nlopt_result.h"

#include "voronoi/new/diagram.h"

#include <nlopt.hpp>
#include <HLBFGS/HLBFGS.h>

#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include<Eigen/IterativeLinearSolvers>

namespace avro
{

namespace voronoi
{

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
  printf("-------------------------------------------------\n");
  printf("%8s|%12s|%12s|%12s\n","iter  "," energy "," gradient "," time (s)");
  printf("-------------------------------------------------\n");
  iteration_ = 0;
  sub_iteration_ = 0;
}

real_t
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

  gnorm = 1.0;
	if      (mode_ == 0) gnorm = calculate_norm(de_dx_.data() , de_dx_.size() );
	else if (mode_ == 1) gnorm = calculate_norm(de_dw_.data() , de_dw_.size() );

	real_t min_vol = *std::min_element( cell_volume_.begin() , cell_volume_.end() );

  printf("%5lu-%3u|%12.3e|%12.3e|%12.3e|%12.3e\n",
    iteration_,sub_iteration_,
    energy_,
    gnorm,
    time_voronoi_,
		min_vol
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

class NewtonStep {

public:
	NewtonStep( const PowerDiagram& diagram ) :
		diagram_(diagram)
	{}

	void build() {

		index_t n = diagram_.sites().nb();
		coord_t dim = diagram_.ambient_dimension();

		// allocate the hessian matrix
		hessian_.setZero();
		hessian_.resize(n,n);
		hessian_.reserve(n*20 ); // TODO count maximum number of facets to estimate size

		// loop through the power cells
		printf("--> building hessian..\n");
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
				avro_assert_msg( pi == k , "pi = %d, k = %lu" , pi , k );
				avro_assert_msg( pi != pj , "pi = %d, pj = %d" , pi , pj );

				const real_t* xj = diagram_.sites()[pj];
				real_t Hij = 0.5*facet.Aij / numerics::distance(xi,xj,dim);

				// set the off-diagonal term
				hessian_.coeffRef(k,pj) = Hij;

				// add the contribution to the diagonal term
				Hii -= Hij;
			}

			// set the diagonal term
			hessian_.coeffRef(k,k) = Hii;
		}
		printf("done\n");

		//std::cout << Eigen::MatrixXd(hessian_) << std::endl;// avro_implement;
	}

	void solve( const std::vector<real_t>& grad , std::vector<real_t>& pk ) {

		index_t n = diagram_.sites().nb();
		Eigen::VectorXd b(n);
		for (index_t k = 0; k < grad.size(); k++) {
			b(k) = -grad[k];
		}

		printf("--> solving linear system..\n");
    //hessian_.makeCompressed();
		//Eigen::SparseLU<Eigen::SparseMatrix<real_t> > solver;
		//igen::BiCGSTAB<Eigen::SparseMatrix<real_t> , Eigen::DiagonalPreconditioner<real_t> > solver; //Eigen::IncompleteLUT<real_t> > solver;
		//Eigen::ConjugateGradient<Eigen::SparseMatrix<real_t> , Eigen::Lower | Eigen::Upper , Eigen::DiagonalPreconditioner<real_t> > solver;
    //Eigen::SparseQR<Eigen::SparseMatrix<real_t> , Eigen::COLAMDOrdering<int> > solver;
    Eigen::LeastSquaresConjugateGradient< Eigen::SparseMatrix<real_t> , Eigen::DiagonalPreconditioner<real_t> > solver;

    solver.setTolerance(1e-8);
    solver.setMaxIterations(10000);

		solver.compute(hessian_);
		if (solver.info() != Eigen::Success) {
			printf("decomposition failed\n");
      avro_assert_not_reached;
		}
		Eigen::VectorXd dw(n);
		dw = solver.solve(b);
		if (solver.info() != Eigen::Success) {
			//std::cout << Eigen::MatrixXd(hessian_) << std::endl;// avro_implement;
			printf("solve failed\n");
      std::cout << eigen_error_message(solver.info()) << std::endl;
      printf("iterations = %d, error = %g",solver.iterations(),solver.error());
      avro_assert_not_reached;
		}
		printf("done\n");

	//	printf("determinant = %g\n",solver.determinant());

		//std::cout << b << std::endl;
		//std::cout << dw << std::endl;

		for (index_t k = 0; k < pk.size(); k++)
			pk[k] = dw(k);
	}

private:
	const PowerDiagram& diagram_;
	Eigen::SparseMatrix<real_t> hessian_;
};

void
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
  real_t gnorm = -1.0;
  real_t a0 = 0.5*( std::min(vmin,mmin) );
  printf("a0 = %g\n",a0);

	for (index_t iter = 0; iter < nb_iter; iter++) {

		// build the hessian matrix
		NewtonStep newton(*this);
		newton.build();

		// determine the descent direction
		std::vector<real_t> pk( n );
		newton.solve( de_dw_ , pk );

    // determine the step size
    real_t alpha = 1.0;

    // copy the original weights, and calculate the current norm of the energy gradient
    std::vector<real_t> x0(x.begin(),x.end());
    real_t gnorm0 = calculate_norm( de_dw_.data() , n );

    sub_iteration_ = 0;
    while (sub_iteration_ < 25) {

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

      printf("alpha = %g\n",alpha);

      if (alpha < 1e-7) avro_assert_not_reached;
    }

    if (gnorm < 1e-12) break;

		iteration_++;
	}

}

} // voronoi

} // avro
