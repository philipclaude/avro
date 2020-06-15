#include "unit_tester.hpp"

#include "common/json.h"
#include "common/tools.h"

#include "library/ckf.h"

#include "mesh/points.h"

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

UT_TEST_SUITE_END( mesh_points_test_suite )
