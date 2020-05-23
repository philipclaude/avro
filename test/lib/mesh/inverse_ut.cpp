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
    UT_CATCH_EXCEPTION( topology.inverse().shell( 0 , topology.points().nb()-1 , S ) );
  }
}
UT_TEST_CASE_END(ckf_nd)

UT_TEST_SUITE_END(mesh_inverse_suite)
