#include "unit_tester.hpp"

#include "graphics/plot.h"
#include "graphics/plotter.h"
#include "graphics/window.h"

#include "mesh/topology.h"

#include "numerics/field.h"

using namespace ursa;
using namespace ursa::graphics;

UT_TEST_SUITE( PlotterSuite )

UT_TEST_CASE( test1 )
{
  Plotter plotter;

  DummyTopology topology;
  Fields fields;
  Window::Plot_ptr plot = std::make_shared<Plot>(topology,&fields);

  plotter.window("main").attach(plot); // make generic later

  plotter.run();
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( PlotterSuite )
