#include "unit_tester.hpp"

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
  //library::objFile topology( "/Users/pcaplan/Desktop/spot_triangulated.obj" );
  library::objFile topology( "/Users/pcaplan/Desktop/suzanne.obj" );

  std::vector<index_t> edges;
  topology.getEdges(edges);

  Topology<Simplex<Lagrange>> topology_edges( topology.vertices() , 1 );
  for (index_t k=0;k<edges.size()/2;k++)
    topology_edges.add( edges.data()+2*k , 2 );

  Window::Plot_ptr plot1 = std::make_shared<Plot>(topology,&window);
  plotter.window("main").attach(plot1);

  plotter.run();
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( PlotterSuite )
