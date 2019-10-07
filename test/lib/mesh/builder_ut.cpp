#include "unit_tester.hpp"

#include "library/samples.h"

#include "mesh/builder.h"

using namespace ursa;

UT_TEST_SUITE( TopologySuite )

UT_TEST_CASE( simplex_tests )
{
  library::TwoTriangles topology;
  Vertices vertices( topology.vertices().dim() );

  Topology<Simplex<Lagrange>> topology_curved( vertices , topology , 3 );

  vertices.print();
  //builder.transfer( topology_curved );

}
UT_TEST_CASE_END( simplex_tests )

UT_TEST_SUITE_END( TopologySuite )
