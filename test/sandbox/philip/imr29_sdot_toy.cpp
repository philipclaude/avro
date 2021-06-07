#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/plots.h"

#include "mesh/facets.h"

#include "numerics/geometry.h"

#include "voronoi/delaunay.h"
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
  index_t nb_points = 1e3;
  index_t nb_qt_iter = 100;
  index_t nb_ot_iter = 200;
  bool lloyd = false;

  coord_t dim = number+1;
  CubeDomain<type> domain(number,dim,2);

  // uniform density
  delaunay::DensityMeasure_Uniform density1(1.0);
  DensityMeasure_Sphere density2(number);
  DensityMeasure_Cone density3(number);

  std::vector<delaunay::DensityMeasure*> measures = { &density3 };//&density2 , &density1 };
  std::vector<std::string> names = {"cone"};//"sphere","uniform"};


  for (index_t imeasure = 0; imeasure < measures.size(); imeasure++)
  {

    delaunay::SemiDiscreteOptimalTransport<type> transport(domain,measures[imeasure]);

    std::string prefix = "/home/pcaplan/Dropbox/research/publications/imr-2021-xxxx/optimal_transport/sdot-" + names[imeasure] + "-dim-" + std::to_string(number) + "-n-" + std::to_string(nb_points);

    printf("prefix = %s\n",prefix.c_str());

    srand(1);
    transport.sample( nb_points );

    transport.weight_max() = 1e-1;
    transport.quad_order() = 3;

    if (names[imeasure] == "uniform") transport.quad_order() = 2;

    transport.save_every( 1e6 , "tmp/sdot-tmp" );

    // first optimize the points
    if (lloyd)
      transport.optimize_points_lloyd(nb_qt_iter);
    else
      transport.optimize_points(nb_qt_iter);

    transport.save_every( 10 , prefix );

    // now optimize the masses
    const std::vector<real_t>& mass = transport.mass();
    real_t mass_total = 0.0;
    for (index_t k = 0; k < mass.size(); k++)
      mass_total += mass[k];
    real_t mass_min = * std::min_element( mass.begin() , mass.end() );
    real_t mass_max = * std::max_element( mass.begin() , mass.end() );
    printf("total mass = %g, min = %g, max = %g, average = %g\n",mass_total,mass_min,mass_max,mass_total/real_t(nb_points));
    std::vector<real_t> nu( nb_points , mass_total / real_t(nb_points) );
    transport.set_nu( nu );
    transport.optimize_weights(nb_ot_iter);

    transport.properties().save( prefix + "-properties.json" );
  }
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( sandbox_semidiscrete_ot_toy )
