#include "numerics/nlopt_result.h"

#include "voronoi/new/diagram.h"

#include <nlopt.hpp>
#include <HLBFGS/HLBFGS.h>

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

  compute();

  if (mode_ == 0) {
    for (index_t i = 0; i < de_dx_.size(); i++)
      grad[i] = de_dx_[i];
  }

  gnorm = calculate_norm(de_dx_.data() , de_dx_.size() );
  printf("%8lu|%12.3e|%12.3e|%12.3e\n",
    iteration_,
    energy_,
    gnorm,
    time_voronoi_
  );

  return energy_;
}

void
PowerDiagram::optimize_points( index_t nb_iter ) {

  mode_ = 0;
  coord_t dim = domain_.number();
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
  #if 1 // this should not be used for general meshes
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
  catch (...) {
    printf("nlopt failed\n");
  }

}

} // voronoi

} // avro
