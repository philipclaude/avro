#include "unit_tester.hpp"

#include "library/ckf.h"

#include "voronoi/delaunay.h"
#include "voronoi/voronoi.h"

UT_TEST_SUITE( voronoi_test_suite )

UT_TEST_CASE( test1 )
{
  coord_t number = 4;
  std::vector<index_t> dims(number,3);
  CKF_Triangulation topology( dims );

  Delaunay delaunay(topology.points());
  delaunay::RestrictedVoronoiDiagram rvd(topology,delaunay);

  rvd.parallel() = true;
  rvd.compute(false);

  //rvd.print();

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( voronoi_test_suite )
