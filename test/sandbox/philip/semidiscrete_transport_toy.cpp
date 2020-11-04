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
  coord_t number = 3;
  index_t nb_points = 1e4;

  coord_t dim = 5;
  CubeDomain<type> domain(number,dim,3);

  delaunay::DensityMeasure_Uniform density(1.0);
  DensityMeasure_Test density2;

  delaunay::SemiDiscreteOptimalTransport<type> transport(domain,&density);
  transport.sample( nb_points );
  transport.optimize_points_lloyd(10);

  std::vector<real_t> m = transport.mass();
  real_t m_total = 0.0;
  for (index_t j = 0; j < m.size() ; j++)
    m_total += m[j];
  real_t m_min = * std::min_element( m.begin() , m.end() );
  real_t m_max = * std::max_element( m.begin() , m.end() );
  printf("mass properties: min = %g, max = %g, total = %g\n",m_min,m_max,m_total);

  // set the mass to the current mass
  #if 0
  std::vector<real_t> mass( nb_points , 1e-4 );
  #else
  std::vector<real_t>& mass = transport.mass();
  mass[0] *= 10;
  #endif
  transport.set_nu( mass );

  // optimize the weights
  transport.set_density( &density );
  transport.optimize_weights(25);

  delaunay::IntegrationSimplices& triangulation = transport.simplices();
  std::shared_ptr<delaunay::TriangulationCells> tc = std::make_shared<delaunay::TriangulationCells>(triangulation);
  triangulation.fields().make("c",tc);

  std::shared_ptr<delaunay::TriangulationElements> te = std::make_shared<delaunay::TriangulationElements>(triangulation);
  triangulation.fields().make("e",te);

  if (number > 3 || (nb_points >= 1e6)) return;
  graphics::Visualizer vis;
  library::Plot<Simplex> point_plot(transport.delaunay());

  vis.add_topology(triangulation);
  vis.add_topology(point_plot);
  //vis.add_topology(transport.diagram());
  vis.run();
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( sandbox_semidiscrete_ot_toy )
