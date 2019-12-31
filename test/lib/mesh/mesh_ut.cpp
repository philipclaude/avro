#include "unit_tester.hpp"

#include "common/tools.h"

#include "mesh/mesh.h"

using namespace avro;

UT_TEST_SUITE( mesh_test_suite )

UT_TEST_CASE( test1 )
{
  Mesh mesh(3);

  std::shared_ptr<Topology<Simplex>> topology = std::make_shared<Topology<Simplex>>(mesh.points(),2);

  mesh.add( topology );

  Topology<Simplex> ts = mesh.template retrieve<Simplex>(0);

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( mesh_test_suite )
