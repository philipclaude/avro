#include "unit_tester.hpp"

#include "common/json.h"
#include "common/tools.h"

#include "library/ckf.h"

#include "mesh/points.h"

#include "numerics/geometry.h"

#include "voronoi/delaunay.h"
#include "voronoi/voronoi.h"

#include <json/json.hpp>

using namespace avro;

UT_TEST_SUITE( mesh_points_test_suite )

UT_TEST_CASE( test1 )
{
  Points points(3);

  coord_t number = 3;
  index_t N = 2;
  std::vector<index_t> dims(number,N);
  CKF_Triangulation topology( dims );

  Delaunay delaunay(topology.points().dim());
  topology.points().copy(delaunay);

  #if 0
  for (index_t k=0;k<delaunay.nb();k++)
  for (coord_t d=0;d<number;d++)
    delaunay[k][d] += 0.5;
  #endif

  printf("running rvd test for %u-simplex mesh with %lu elements and %lu delaunay vertices\n",number,topology.nb(),delaunay.nb());

  delaunay::RestrictedVoronoiDiagram rvd(topology,delaunay);
  rvd.parallel() = true;

  rvd.compute(true);

  std::vector<index_t> idx1;
  rvd.points().duplicates(idx1);

  std::vector<index_t> idx2;
  rvd.points().duplicates(idx2,rvd.points().incidence() );

  UT_ASSERT_EQUALS( idx1.size() , idx2.size() );

  for (index_t k=0;k<idx1.size();k++)
  {
    //UT_ASSERT_EQUALS( idx1[k] , idx2[k] );
  }

  json J;
  rvd.points().to_json( J );

  rvd.points().print(0,true);

  rvd.points().dump( "tmp/points.dat" );

}
UT_TEST_CASE_END( test1 )

UT_TEST_CASE( points_move )
{
  coord_t number = 3;
  index_t N = 4;
  std::vector<index_t> dims(number,N);
  CKF_Triangulation topology( dims );

  Points& points = topology.points();
  Points points0 = points;
  points.print();

  index_t k1 = 0;
  for (index_t k=1;k<points0.nb();k++)
  {
    points.print(k);
    points.move_to( k , k1 );
    points.print(k1);
    points0.print(k);

    real_t d = numerics::distance( points0[k] , points[k1] , points.dim() );
    UT_ASSERT_NEAR( d , 0.0 , 1e-12 );

    // make sure all the other vertices are the same
    for (index_t j=k+1;j<points0.nb();j++)
    {
      d = numerics::distance( points0[j] , points[j] , points.dim() );
      UT_ASSERT_NEAR( d , 0.0 , 1e-12 );
    }

    k1++;
  }

  points.print();


}
UT_TEST_CASE_END( points_move )

UT_TEST_SUITE_END( mesh_points_test_suite )
