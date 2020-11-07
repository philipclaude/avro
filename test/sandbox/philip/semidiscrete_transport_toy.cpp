#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/plots.h"

#include "mesh/facets.h"

#include "numerics/geometry.h"

#include "voronoi/delaunay.h"
#include "voronoi/optimal_transport.h"

#include <fstream>

UT_TEST_SUITE( sandbox_semidiscrete_ot_toy )

class DensityMeasure_Test : public delaunay::DensityMeasure
{
public:
  real_t evaluate( index_t elem , const real_t* xref , const real_t* x ) const
  {
    return 1e1*( 1 + sin(2*M_PI*x[0])*sin(2*M_PI*x[1]) );
  }
};

UT_TEST_CASE( test1 )
{
  typedef Polytope type;
  //typedef Simplex type;
  coord_t number = 2;
  index_t nb_points = 1e2;

  coord_t dim = 3;
  CubeDomain<type> domain(number,dim,2);

  delaunay::DensityMeasure_Uniform density(1.0);
  DensityMeasure_Test density2;

  delaunay::SemiDiscreteOptimalTransport<type> transport(domain,&density);
  transport.sample( nb_points );

  #if 0
  transport.optimize_points(30);

  std::vector<real_t> m = transport.mass();
  real_t m_total = 0.0;
  for (index_t j = 0; j < m.size() ; j++)
    m_total += m[j];
  real_t m_min = * std::min_element( m.begin() , m.end() );
  real_t m_max = * std::max_element( m.begin() , m.end() );
  printf("mass properties: min = %g, max = %g, total = %g\n",m_min,m_max,m_total);

  // set the mass to the current mass
  #if 1
  std::vector<real_t> mass( nb_points , 1.0/nb_points );
  #else
  std::vector<real_t>& mass = transport.mass();
  mass[0] *= 10;
  #endif
  transport.set_nu( mass );

  // optimize the weights
  //transport.set_density( &density2 );
  //transport.optimize_weights(30);
  #else
  transport.generate_bluenoise();

  json J;
  std::vector<real_t> radius( transport.delaunay().nb() );
  std::vector<real_t> x( transport.delaunay().nb()*number );
  std::vector<real_t> center(number,0.5);
  index_t i = 0;
  for (index_t k = 0; k < transport.delaunay().nb(); k++)
  {
    radius[k] = numerics::distance( center.data() , transport.delaunay()[k] , number );
    for (coord_t d = 0; d < number; d++)
      x[i++] = transport.delaunay()[k][d];
  }
  J["r"] = radius;
  J["x"] = x;

  std::ofstream output("bluenoise-dim"+std::to_string(number)+"-n"+std::to_string(transport.delaunay().nb())+".json");
  output << std::setw(4) << J << std::endl;

  #endif

  if (number > 3 || (nb_points >= 1e6)) return;
  delaunay::IntegrationSimplices& triangulation = transport.simplices();
  std::shared_ptr<delaunay::TriangulationCells> tc = std::make_shared<delaunay::TriangulationCells>(triangulation);
  triangulation.fields().make("c",tc);

  std::shared_ptr<delaunay::TriangulationElements> te = std::make_shared<delaunay::TriangulationElements>(triangulation);
  triangulation.fields().make("e",te);

  graphics::Visualizer vis;
  library::Plot<Simplex> point_plot(transport.delaunay());

  vis.add_topology(triangulation);
  vis.add_topology(point_plot);
  //vis.add_topology(transport.diagram());
  vis.run();
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( sandbox_semidiscrete_ot_toy )
