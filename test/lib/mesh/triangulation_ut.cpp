#include "unit_tester.hpp"

#include "graphics/plot.h"
#include "graphics/plotter.h"
#include "graphics/window.h"

#include "library/ckf.h"

#include "mesh/triangulation.h"

#include "voronoi/delaunay.h"
#include "voronoi/voronoi.h"

using namespace avro;
using namespace avro::graphics;

UT_TEST_SUITE( triangulation_suite )

UT_TEST_CASE( simplex_topologies )
{

}
UT_TEST_CASE_END( simplex_topologies )

UT_TEST_CASE( voronoi_tests )
{
  CKF_Triangulation topology( {2,2} );

  Delaunay z(topology.points());

  delaunay::RestrictedVoronoiDiagram rvd(topology,z);
  rvd.compute();

  Triangulation<Polytope> triangulation(rvd);
  triangulation.extract();

  triangulation.print();
  triangulation.points().print();

  Plotter plotter;
  Window& window = plotter.window("main");
  Window::Plot_ptr plot1 = std::make_shared<Plot>(triangulation,&window);
  plotter.window("main").attach(plot1);
  plotter.run();

}
UT_TEST_CASE_END( voronoi_tests )

UT_TEST_SUITE_END( triangulation_suite )
