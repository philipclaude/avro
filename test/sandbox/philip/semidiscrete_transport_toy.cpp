#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"

#include "mesh/facets.h"

#include "numerics/geometry.h"

#include "voronoi/delaunay.h"
#include "voronoi/optimal_transport.h"

#include <fstream>

UT_TEST_SUITE( sandbox_semidiscrete_ot_toy )

UT_TEST_CASE( test1 )
{
  typedef Polytope type;
  //typedef Simplex type;
  coord_t number = 2;
  index_t nb_points = 1e4;

  coord_t dim = number;
  CubeDomain<type> domain(dim,10);

  delaunay::DensityMeasure_Example density;

  delaunay::SemiDiscreteOptimalTransport<type> transport(domain,density);
  transport.sample( nb_points );

  transport.evaluate();

  delaunay::IntegrationSimplices& triangulation = transport.simplices();
  std::shared_ptr<delaunay::TriangulationCells> tc = std::make_shared<delaunay::TriangulationCells>(triangulation);
  triangulation.fields().make("c",tc);

  std::shared_ptr<delaunay::TriangulationElements> te = std::make_shared<delaunay::TriangulationElements>(triangulation);
  triangulation.fields().make("e",te);

  if (number > 3 || (nb_points >= 1e6)) return;
  graphics::Visualizer vis;
  vis.add_topology(triangulation);
  //vis.add_topology(diagram);
  vis.run();
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( sandbox_semidiscrete_ot_toy )
