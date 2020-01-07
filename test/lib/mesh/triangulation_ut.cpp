#include "unit_tester.hpp"

#include "graphics/plot.h"
#include "graphics/plotter.h"
#include "graphics/window.h"

#include "library/ckf.h"
#include "library/samples.h"

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
  CKF_Triangulation topology( {4,4,4} );
  library::RegularPolygon polygon(6);

  Delaunay z(topology.points());
  for (index_t k=0;k<z.nb();k++)
  for (index_t d=0;d<z.dim();d++)
    z[k][d] = random_within( 0.0 , 1.0 );

  delaunay::RestrictedVoronoiDiagram rvd(topology,z);
  rvd.compute();

  Plotter plotter;
  Window& window = plotter.window("main");
  Window::Plot_ptr plot1 = std::make_shared<Plot>(rvd,&window);
  Window::Plot_ptr plot2 = std::make_shared<Plot>(topology,&window);
  Window::Plot_ptr plot3 = std::make_shared<Plot>(polygon,&window);
  plotter.window("main").attach(plot1);
  plotter.run();
}
UT_TEST_CASE_END( voronoi_tests )

UT_TEST_SUITE_END( triangulation_suite )
