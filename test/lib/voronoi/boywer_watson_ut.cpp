#include "unit_tester.hpp"

#include "common/tools.h"

#include "graphics/application.h"
#include "graphics/user_interface.h"
#include "graphics/window.h"

#include "library/ckf.h"

#include "mesh/points.h"

#include "voronoi/algorithms.h"

using namespace avro;

UT_TEST_SUITE( voronoi_boywer_watson_suite )

UT_TEST_CASE( test0)
{

  coord_t dim = 2;

  for (index_t i=3;i<4;i++)
  {
    index_t nb_points = pow(10,i);

    Points points(dim);
    for (index_t k=0;k<nb_points;k++)
    {
      std::vector<real_t> x(dim);
      for (coord_t d=0;d<dim;d++)
        x[d] = random_within( -1.0 , 1.0 );
      points.create(x.data());
    }

    BowerWatson delaunay(points);
    delaunay.compute();

    #if 1
    graphics::Visualizer vis;
    vis.add_topology(delaunay);

    std::shared_ptr<graphics::Widget> toolbar = std::make_shared<graphics::Toolbar>(vis.main_window(),vis);
    vis.main_window().interface().add_widget( toolbar );

    vis.run();
    #endif

  }

}
UT_TEST_CASE_END( test0 )


UT_TEST_CASE( test1 )
{

  coord_t dim = 2;

  //FILE *fid = fopen("delaunay_timing.dat","w");

  for (index_t i=0;i<3;i++)
  {
    index_t nb_points = pow(10,i);

    Points points(dim);
    for (index_t k=0;k<nb_points;k++)
    {
      std::vector<real_t> x(dim);
      for (coord_t d=0;d<dim;d++)
        x[d] = random_within( -1.0 , 1.0 );
      points.create(x.data());
    }

    BowerWatson delaunay(points);

    clock_t t0 = clock();
    delaunay.compute();
    clock_t t1 = clock();
    UNUSED(t0);
    UNUSED(t1);

    //fprintf(fid,"%lu %g\n",nb_points,real_t(t1-t0)/CLOCKS_PER_SEC);
  }
  //fclose(fid);

  #if 0
  graphics::Visualizer vis;
  vis.add_topology(delaunay);

  std::shared_ptr<graphics::Widget> toolbar = std::make_shared<graphics::Toolbar>(vis.main_window(),vis);
  vis.main_window().interface().add_widget( toolbar );

  vis.run();
  #endif

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( voronoi_boywer_watson_suite )
