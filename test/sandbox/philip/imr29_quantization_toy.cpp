#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/plots.h"

#include "mesh/facets.h"

#include "numerics/geometry.h"

#include "voronoi/optimal_transport.h"

#include "measures.h"
#include "visualize.h"

#include <fstream>
#include <iomanip>

UT_TEST_SUITE( sandbox_semidiscrete_ot_toy )

UT_TEST_CASE( test1 )
{
  typedef Polytope type;
  coord_t number = 4;
  index_t nb_points = 1e4;
  index_t nb_iter = 100;

  coord_t dim = number+1;
  CubeDomain<type> domain(number,dim,2);

  // uniform density
  delaunay::DensityMeasure_Uniform density1(1.0);

  // gaussian
  vecd<real_t> mu(number);
  symd<real_t> sigma(number,number);
  sigma = 0;
  for (coord_t d = 0; d < number; d++)
  {
    mu(d) = 0.5;
    sigma(d,d) = 0.02;
  }
  DensityMeasure_Gaussian density2(mu,sigma);

  // shock
  DensityMeasure_Cone density3(number);


  std::vector<delaunay::DensityMeasure*> measures = { &density3 , &density2 , &density1 };
  std::vector<std::string> names = {"cone","gaussian","uniform"};

  //std::vector<delaunay::DensityMeasure*> measures = { &density2 };
  //std::vector<std::string> names = {"shock"};

  for (index_t ialg = 0; ialg < 2 ; ialg++)
  {

    bool lloyd = false;
    if (ialg == 0) lloyd = false;
    else lloyd = true;

    for (index_t imeasure = 0; imeasure < measures.size(); imeasure++)
    {

      delaunay::SemiDiscreteOptimalTransport<type> transport(domain,measures[imeasure]);

      std::string prefix = "/home/pcaplan/Dropbox/research/publications/imr-2021-xxxx/qntz-" + names[imeasure] + "-" + ((lloyd) ? "lloyd" : "lbfgs") + "-dim-" + std::to_string(number) + "-n-" + std::to_string(nb_points);

      printf("prefix = %s\n",prefix.c_str());

      //if (names[imeasure] == "shock") nb_points = 1e4;
      //else nb_points = 1e3;

      srand(1);
      transport.sample( nb_points );

      transport.save_every( 10 , prefix );

      transport.weight_max() = 1e-1;
      transport.quad_order() = 3;

      if (names[imeasure] == "uniform") transport.quad_order() = 2;

      if (lloyd)
        transport.optimize_points_lloyd(nb_iter);
      else
        transport.optimize_points(nb_iter);

      transport.properties().save( prefix + "-properties.json" );
    }
  }
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( sandbox_semidiscrete_ot_toy )
