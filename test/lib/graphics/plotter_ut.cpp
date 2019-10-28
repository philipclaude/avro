#include "unit_tester.hpp"

#include "graphics/interface.h"
#include "graphics/plot.h"
#include "graphics/plotter.h"
#include "graphics/window.h"

#include "library/obj.h"
#include "library/samples.h"

#include "mesh/topology.h"

using namespace ursa;
using namespace ursa::graphics;

UT_TEST_SUITE( PlotterSuite )

UT_TEST_CASE( test1 )
{
  Plotter plotter;

  Window& window = plotter.window("main");

  printf("reading obj...\n");

  //library::TwoTriangles topology;
  library::objFile topology( "/Users/pcaplan/Google Drive/library/models/obj/spot.obj" );

  Window::Plot_ptr plot1 = std::make_shared<Plot>(topology,&window);
  plotter.window("main").attach(plot1);

  BasicInterface basic( window );
  PlotTree tree(window);
  window.set_interface(&tree);

  plotter.run();
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( PlotterSuite )
