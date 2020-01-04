#include "unit_tester.hpp"

#include "library/ckf.h"

#include "voronoi/delaunay.h"
#include "voronoi/voronoi.h"

UT_TEST_SUITE( voronoi_test_suite )

UT_TEST_CASE( test1 )
{
  coord_t number = 3;
  std::vector<index_t> dims(number,7);
  CKF_Triangulation topology( dims );

  Delaunay delaunay(topology.points().dim());
  topology.points().copy(delaunay);
  for (index_t k=0;k<delaunay.nb();k++)
  for (coord_t d=0;d<number;d++)
    delaunay[k][d] += 0.5;

  delaunay::RestrictedVoronoiDiagram rvd(topology,delaunay);

  rvd.parallel() = true;
  rvd.compute(true);

  //rvd.print();

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( voronoi_test_suite )
