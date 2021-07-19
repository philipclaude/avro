//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"

#include "mesh/decomposition.h"

#include "voronoi/delaunay.h"
#include "voronoi/voronoi.h"

using namespace avro;
using namespace avro::graphics;

UT_TEST_SUITE( voronoi_test_suite )

#if 1
UT_TEST_CASE( line_test )
{
  coord_t number = 1;

  index_t nx = 5;
  Points points(2);
  real_t x[2] = { 0. , 0. };
  index_t e[2] = {0,0};
  Topology<Simplex> topology(points,1);
  for (index_t k=0;k<nx;k++)
  {
    x[1] = real_t(k)/real_t(nx-1.0);
    points.create(x);

    if (k<nx-1)
    {
      e[0] = k;
      e[1] = k+1;
      topology.add( e , 2 );
    }
  }
  points.print(true);

  Delaunay delaunay(points.dim());
  #if 0
  index_t np = 4;
  x[0] = 0.0;
  delaunay.create(x);
  for (index_t k=0;k<np;k++)
  {
    x[0] = random_within(0.1,0.9);
    delaunay.create(x);
  }
  x[0] = 1.0;
  delaunay.create(x);
  #else
  topology.points().copy(delaunay);
  #endif
  delaunay.print(true);

  printf("running rvd test for %u-simplex mesh with %lu elements and %lu delaunay vertices\n",number,topology.nb(),delaunay.nb());

  delaunay::RestrictedVoronoiDiagram rvd(topology,delaunay);
  rvd.parallel() = true;

  rvd.compute(true);

  printf("rvd points:\n");
  rvd.points().print();

  printf("rvd:\n");
  rvd.Table<index_t>::print();
  printf("vfm:\n");
  rvd.points().incidence().print();

  Topology<Simplex> dt( delaunay , number );
  printf("extracting dt:\n");
  rvd.extract(dt);
  printf("delaunay triangulation has %lu simplices\n",dt.nb());

  Viewer vis;

  //vis.add(topology);
  vis.add(rvd);
  vis.add(dt);
  //vis.add(T);

  //vis.run(AVRO_FULL_UNIT_TEST);
}
UT_TEST_CASE_END( line_test )
#endif

UT_TEST_CASE( test0 )
{
  coord_t number = 2;
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

  Topology<Simplex> dt( delaunay , number );
  rvd.extract(dt);
  printf("delaunay triangulation has %lu simplices\n",dt.nb());

  Viewer vis;

  //vis.add(topology);
  vis.add(rvd);
  //vis.add(T);

  //vis.run(AVRO_FULL_UNIT_TEST);
}
UT_TEST_CASE_END( test0 )

UT_TEST_CASE( test1 )
{
  for (coord_t number=2;number<=4;number++)
  {
    for (index_t N=2;N<=4;N++)
    {

      if (number==4 && N >= 4) break; // very slow!

      std::vector<index_t> dims(number,N);
      CKF_Triangulation topology( dims );

      Delaunay delaunay(topology.points().dim());
      topology.points().copy(delaunay);

      printf("running rvd test for %u-simplex mesh with %lu elements and %lu delaunay vertices\n",number,topology.nb(),delaunay.nb());

      delaunay::RestrictedVoronoiDiagram rvd(topology,delaunay);
      rvd.parallel() = true;

      // test 1: sites at mesh points (no need for exact precision)
      rvd.compute(false);

      // test 2: sites offset to test exact precision
      for (index_t k=0;k<delaunay.nb();k++)
      for (coord_t d=0;d<number;d++)
        delaunay[k][d] += 0.5;

      // must run with exact precision
      rvd.compute(true);

      SimplicialDecomposition<Polytope> simplices(rvd);
      simplices.extract();
      std::vector<index_t> S,P;
      simplices.get_simplices( rvd.number(),S,P );

      Topology<Simplex> T(simplices.points(),rvd.number());
      for (index_t k=0;k<S.size()/(rvd.number()+1);k++)
        T.add( &S[k*(rvd.number()+1)],rvd.number()+1);
      printf("--> volume = %g\n",T.volume());

      Topology<Simplex> dt( delaunay , number );
      rvd.extract(dt);
      printf("delaunay triangulation has %lu simplices\n",dt.nb());
    }
  }

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( voronoi_test_suite )
