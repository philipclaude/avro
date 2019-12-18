#include "adaptation/adapt.h"

#include <set>

namespace luna
{

template<typename type>
void
AdaptThread<type>::smooth( MetricField<type>& metric , const index_t niter , FILE* fid )
{
  std::vector<index_t> elems;
  std::vector<index_t> N;

  smoother_.nb_parameter_tests() = 0;
  smoother_.nb_parameter_rejections() = 0;
  smoother_.nb_geometry() = 0;

  topology_.evaluate(metric);
  real Q0 = topology_.worst();

  // loop over smoothing iterations
  printf("-> performing vertex smoothing:\n");
  for (index_t iter=0;iter<niter;iter++)
  {

    smoother_.delta() = 0.0;
    smoother_.delta_min() = 1e20;
    smoother_.delta_max() = -1.;

    smoother_.resetRejections();
    smoother_.nb_accepted() = 0;

    smoother_.objective() = 0.0;

    // loop through the vertices
    for (index_t k=0;k<topology_.vertices().nb();k++)
    {
      if (k<topology_.vertices().nb_ghost()) continue;
      smoother_.apply( k , metric , Q0 );
    }
    smoother_.objective() /= (topology_.vertices().nb()/topology_.vertices().dim());

    printf("\titer[%lu]: dx = %3.2e, min = %3.2e, max = %3.2e -> accepted %lu\n",iter,smoother_.delta()/topology_.vertices().nb(),smoother_.delta_min(),smoother_.delta_max(),
          smoother_.nb_accepted());
    printf("\t\tnb_visibility_rej = %lu, nb_implied_metric_rej = %lu, nb_enlarged_rej = %lu/%lu, Navg = %lu\n",
      smoother_.nb_visibility_rejections(),smoother_.nb_implied_metric_rejections(),
      smoother_.nb_enlarged_rejections(),smoother_.nb_geometry(),index_t(smoother_.Ntot()/topology_.vertices().nb()));
      printf("\t\tobjective = %1.12e, nb_error = %lu, nb_interp_outside = %lu\n",smoother_.objective(),smoother_.nb_zero_valency(),smoother_.nb_interpolated_outside());

    if (fid!=NULL)
      fprintf(fid,"%lu %1.12e\n",iter,smoother_.objective());
  }
  printf("\tdone %lu iterations of smoothing.\n\tnb_elem = %lu, nb_vert = %lu. parameter rej = (%lu/%lu)\n",
    niter,topology_.nb(),topology_.vertices().nb(),
    smoother_.nb_parameter_rejections(),smoother_.nb_parameter_tests());
}

} // luna
