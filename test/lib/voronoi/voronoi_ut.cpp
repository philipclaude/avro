#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"

#include "mesh/decomposition.h"

#include "voronoi/delaunay.h"
#include "voronoi/voronoi.h"

using namespace avro;
using namespace avro::graphics;

UT_TEST_SUITE( voronoi_test_suite )

UT_TEST_CASE( test0 )
{
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

  rvd.Table<index_t>::print();
  rvd.points().incidence().print();

  SimplicialDecomposition<Polytope> simplices(rvd);
  simplices.extract();
  std::vector<index_t> S,P;
  simplices.get_simplices( rvd.number(),S,P );

  Topology<Simplex> T(simplices.points(),rvd.number());
  for (index_t k=0;k<S.size()/(rvd.number()+1);k++)
    T.add( &S[k*(rvd.number()+1)],rvd.number()+1);
  printf("--> volume = %g\n",T.volume());

  //T.Table<index_t>::print();

  Visualizer vis;

  //vis.add_topology(topology);
  vis.add_topology(rvd);
  //vis.add_topology(T);

  vis.run();
}
UT_TEST_CASE_END( test0 )

UT_TEST_CASE( test1 )
{
  return;
  for (coord_t number=2;number<=4;number++)
  {
    for (index_t N=2;N<=4;N++)
    {

      if (number==4 && N > 3) break; // very slow!

      std::vector<index_t> dims(number,N);
      CKF_Triangulation topology( dims );

      Delaunay delaunay(topology.points().dim());
      topology.points().copy(delaunay);

      printf("running rvd test for %u-simplex mesh with %lu elements and %lu delaunay vertices\n",number,topology.nb(),delaunay.nb());

      delaunay::RestrictedVoronoiDiagram rvd(topology,delaunay);
      rvd.parallel() = true;

      // test 1: sites at mesh points
      rvd.compute(true);
      rvd.compute(false);

      // test 2: sites offset to test exact precision
      for (index_t k=0;k<delaunay.nb();k++)
      for (coord_t d=0;d<number;d++)
        delaunay[k][d] += 0.5;

      rvd.compute(true);
      rvd.compute(false);

      SimplicialDecomposition<Polytope> simplices(rvd);
      simplices.extract();
      std::vector<index_t> S,P;
      simplices.get_simplices( rvd.number(),S,P );

      Topology<Simplex> T(simplices.points(),rvd.number());
      for (index_t k=0;k<S.size()/(rvd.number()+1);k++)
        T.add( &S[k*(rvd.number()+1)],rvd.number()+1);
      printf("--> volume = %g\n",T.volume());
    }
  }

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( voronoi_test_suite )
