#include "unit_tester.hpp"

#include "library/ckf.h"

#include "mesh/partition.h"

using namespace avro;

UT_TEST_SUITE( mesh_partition_test_suite )

UT_TEST_CASE( test1 )
{
  coord_t number = 3;

  std::vector<index_t> dims(number,10);
  CKF_Triangulation topology( dims );
  //topology.close();
  topology.neighbours().compute();

  Partition<Simplex> partition(topology);
  partition.compute(4);
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( mesh_partition_test_suite )
