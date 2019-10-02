#include "unit_tester.hpp"

#include "graphics/plot.h"
#include "graphics/plotter.h"
#include "graphics/window.h"

#include "library/samples.h"

#include "mesh/topology.h"

#include "numerics/field.h"

using namespace ursa;
using namespace ursa::graphics;

UT_TEST_SUITE( PlotterSuite )

UT_TEST_CASE( test1 )
{
  Plotter plotter;

  Window& window = plotter.window("main");

  library::TwoTriangles topology;
  Window::Plot_ptr plot1 = std::make_shared<Plot>(topology,&window);
  plotter.window("main").attach(plot1);

  Window::Plot_ptr plot2 = std::make_shared<Plot>(topology.edges(),&window);
  //plotter.window("main").attach(plot2);

  plotter.run();
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( PlotterSuite )
