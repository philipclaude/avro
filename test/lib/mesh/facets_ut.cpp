#include "unit_tester.hpp"

#include "library/ckf.h"

#include "mesh/facets.h"
#include "mesh/topology.h"

using namespace avro;

UT_TEST_SUITE(facet_test_suite)

UT_TEST_CASE(ckf_nd)
{
  for (coord_t dim=4;dim<=4;dim++)
  {
    index_t N = 6;
    std::vector<index_t> dims(dim,N);
    CKF_Triangulation topology( dims );

    Facets facets( topology );

    // with a boundary, some facets will not get counted
    clock_t t0 = clock();
    facets.compute();
    clock_t t1 = clock();

    printf("facet time = %g seconds\n",real_t(t1-t0)/CLOCKS_PER_SEC);
    UT_ASSERT_EQUALS( facets.check() , false );

    // close the topology, meaning the check should not fail
    topology.close();
    facets.clear();
    facets.compute();
    UT_ASSERT_EQUALS( facets.check() , true );
  }
}
UT_TEST_CASE_END(ckf_nd)

UT_TEST_SUITE_END(facet_test_suite)
