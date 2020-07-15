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

#include "library/ckf.h"

#include "mesh/decomposition.h"
#include "mesh/inverse.h"
#include "mesh/topology.h"

using namespace avro;

UT_TEST_SUITE(mesh_inverse_suite)

bool
contains( index_t value , const std::vector<index_t>& x )
{
  for (index_t j=0;j<x.size();j++)
    if (x[j]==value) return true;
  return false;
}


UT_TEST_CASE(ckf_nd)
{
  for (coord_t number=2;number<=4;number++)
  {

    std::vector<index_t> dims(number,4);

    CKF_Triangulation topology( dims );

    topology.orient();
    topology.close();
    topology.neighbours().compute();
    topology.inverse().build();

    topology.inverse().print();

    UT_ASSERT( topology.inverse().check() );

    SimplicialDecomposition<Simplex> decomposition(topology);
    decomposition.extract();

    for (index_t j=0;j<number;j++)
    {
      std::vector<index_t> simplices;
      std::vector<index_t> parents;

      decomposition.get_simplices( j , simplices , parents );

      index_t nv = j +1;
      index_t ns = simplices.size()/nv;
      for (index_t k=0;k<ns;k++)
      {
        std::vector<index_t> shell;

        if (j==0) topology.inverse().ball( simplices[k] , shell );
        if (j==1) topology.inverse().shell( simplices[2*k] , simplices[2*k+1] , shell );
        if (j==2) topology.inverse().shell( simplices[3*k] , simplices[3*k+1], simplices[3*k+2] , shell );

        bool result = contains( parents[k] , shell );
        UT_ASSERT( result );
      }

    }
    topology.inverse().print();

    printf("checking inverse topology the long way...\n");
    UT_ASSERT( topology.inverse().check() );
    printf("done.\n");

    InverseTopology<Simplex> inverse_copy( topology );
    inverse_copy.copy( topology.inverse() );

    // create a dummy vertex not conected to any elements
    std::vector<real_t> x(number,0);
    topology.points().create( x.data() );
    topology.inverse().create(1);
    std::vector<index_t> S;
    //UT_CATCH_EXCEPTION( topology.inverse().shell( 0 , topology.points().nb()-1 , S ) );
  }
}
UT_TEST_CASE_END(ckf_nd)

UT_TEST_CASE( non_manifold_topology )
{
  coord_t number = 2;
  coord_t dim = number;

  Points points(dim);

  real_t x0[2] = {0,0};
  real_t x1[2] = {1,0};
  real_t x2[2] = {1,1};
  real_t x3[2] = {2,1};
  real_t x4[2] = {1,2};
  real_t x5[2] = {2,2};
  points.create(x0);
  points.create(x1);
  points.create(x2);
  points.create(x3);
  points.create(x4);
  points.create(x5);

  Topology<Simplex> topology(points,number);
  index_t t0[3] = {0,1,2};
  index_t t1[3] = {2,5,4};
  index_t t2[3] = {2,3,5};
  topology.add(t0,3);
  topology.add(t1,3);
  topology.add(t2,3);

  topology.close();
  topology.build_structures();

  std::vector<index_t> C;
  std::vector<index_t> edges;
  topology.get_edges(edges);
  for (index_t k=0;k<edges.size()/2;k++)
  {
    topology.inverse().shell( edges[2*k],edges[2*k+1] , C );
    print_inline(C);
  }

  topology.Table<index_t>::print();
  //topology.neighbours().print();

}
UT_TEST_CASE_END( non_manifold_topology )

UT_TEST_SUITE_END(mesh_inverse_suite)
