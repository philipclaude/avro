#include "unit_tester.hpp"

#include "graphics/interface.h"
#include "graphics/plot.h"
#include "graphics/plotter.h"
#include "graphics/window.h"

#include "library/ckf.h"
#include "library/obj.h"
#include "library/samples.h"

#include "mesh/topology.h"

using namespace luna;
using namespace luna::graphics;

UT_TEST_SUITE( PlotterSuite )

UT_TEST_CASE( test1 )
{
  Plotter plotter;

  Window& window = plotter.window("main");

  //library::TwoTriangles topology;
  //library::objFile topology( "/Users/pcaplan/Google Drive/library/models/obj/spot.obj" );

  #if 0
  CKF_Triangulation topology0( {10,10,10} );
  Points points(3);
  for (index_t k=0;k<topology0.points().nb();k++)
  {
    std::vector<real_t> x(3,0);
    x.assign(topology0.points()[k],topology0.points()[k]+2);
    points.create(x.data());
  }
  Topology<Simplex> topology(points,2);
  for (index_t k=0;k<topology0.nb();k++)
    topology.add( topology0(k) , topology0.nv(k) );
  #else
  CKF_Triangulation topology( {10,10,10} );
  #endif

  Window::Plot_ptr plot1 = std::make_shared<Plot>(topology,&window);
  plotter.window("main").attach(plot1);

  BasicInterface basic( window );
  PlotTree tree(window);
  window.set_interface(&tree);

  plotter.run();
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( PlotterSuite )
