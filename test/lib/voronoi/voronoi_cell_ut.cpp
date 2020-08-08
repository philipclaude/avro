#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"

#include "voronoi/delaunay.h"
#include "voronoi/voronoi_cell.h"

using namespace avro;

UT_TEST_SUITE( voronoi_cell_test_suite )

UT_TEST_CASE( test_2d )
{
  //return;
  GEO::PCK::initialize();

  coord_t number = 2;
  coord_t dim = number;

  Points points(dim);
  Topology<Polytope> domain(points,number);
  std::vector<index_t> square = {0,1,2,3};
  domain.add( square.data() , square.size() );

  real_t x0[2] = {0,0};
  real_t x1[2] = {1,0};
  real_t x2[2] = {1,1};
  real_t x3[2] = {0,1};

  int f0 = -1; // 0-1
  int f1 = -2; // 1-2
  int f2 = -3; // 2-3
  int f3 = -4; // 3-0

  int b0[2] = {f0,f3};
  int b1[2] = {f0,f1};
  int b2[2] = {f1,f2};
  int b3[2] = {f2,f3};

  points.create( x0 );
  points.create( x1 );
  points.create( x2 );
  points.create( x3 );

  points.incidence().add( b0 , 2 );
  points.incidence().add( b1 , 2 );
  points.incidence().add( b2 , 2 );
  points.incidence().add( b3 , 2 );

  // create random delaunay vertices
  Delaunay delaunay( dim );
  #if 1
  index_t nb_points = 1e2;
  std::vector<real_t> x(dim,0.);
  for (index_t k=0;k<nb_points;k++)
  {
    for (index_t d=0;d<dim;d++)
      x[d] = random_within(0.,1.);
    delaunay.create(x.data());
  }
  #else
  std::vector<index_t> dims(number,10);
  CKF_Triangulation ckf(dims);
  ckf.points().copy(delaunay);
  #endif

  //NearestNeighbours neighbours(delaunay,100);
  //neighbours.compute();
  //delaunay::VoronoiCell cell(0,delaunay,neighbours,domain,true);
  //cell.compute();

  delaunay::VoronoiDiagram diagram( delaunay , domain );
  diagram.compute(false);

  graphics::Visualizer vis;
  vis.add_topology(diagram);

  vis.run();
}
UT_TEST_CASE_END( test_2d )

UT_TEST_CASE( test_2d_simplex )
{
//  return;

  GEO::PCK::initialize();

  coord_t number = 3;
  coord_t dim = number;

  CKF_Triangulation domain( {2,2,2} );

  // create random delaunay vertices
  Delaunay delaunay( dim );
  #if 1
  index_t nb_points = 1e2;
  std::vector<real_t> x(dim,0.);
  for (index_t k=0;k<nb_points;k++)
  {
    for (index_t d=0;d<dim;d++)
      x[d] = random_within(0.,1.);
    delaunay.create(x.data());
  }
  #else
  std::vector<index_t> dims(number,10);
  CKF_Triangulation ckf(dims);
  ckf.points().copy(delaunay);
  #endif

  delaunay::VoronoiDiagram diagram( delaunay , domain , true );
  diagram.compute(true);

  graphics::Visualizer vis;
  vis.add_topology(diagram);

  vis.run();
}
UT_TEST_CASE_END( test_2d_simplex )

UT_TEST_SUITE_END( voronoi_cell_test_suite )
