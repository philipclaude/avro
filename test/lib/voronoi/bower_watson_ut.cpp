#include "unit_tester.hpp"

#include "common/tools.h"

#include "graphics/application.h"
#include "graphics/user_interface.h"
#include "graphics/window.h"

#include "library/ckf.h"

#include "mesh/points.h"

#include "voronoi/algorithms.h"

using namespace avro;

UT_TEST_SUITE( voronoi_bower_watson_suite )

UT_TEST_CASE( test1 )
{

  index_t nb_points = 10000;
  coord_t dim = 2;

  Points points(dim);
  for (index_t k=0;k<nb_points;k++)
  {
    std::vector<real_t> x(dim);
    for (coord_t d=0;d<dim;d++)
      x[d] = random_within( -1.0 , 1.0 );
    points.create(x.data());
  }

  #if 0
  points.clear();
  std::vector<index_t> dims(dim,4);
  CKF_Triangulation triangulation(dims);
  triangulation.points().copy(points);
  #endif


  BowerWatson delaunay(points);

  delaunay.compute();

  //points.print();
  //delaunay.points().print();

  graphics::Visualizer vis;
  vis.add_topology(delaunay);

  std::shared_ptr<graphics::Widget> toolbar = std::make_shared<graphics::Toolbar>(vis.main_window(),vis);
  vis.main_window().interface().add_widget( toolbar );

  vis.run();

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( voronoi_bower_watson_suite )
