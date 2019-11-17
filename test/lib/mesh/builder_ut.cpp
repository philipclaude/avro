#include "unit_tester.hpp"

#include "library/samples.h"

#include "mesh/builder.h"

using namespace luna;

UT_TEST_SUITE( TopologySuite )

UT_TEST_CASE( simplex_tests )
{
  library::TwoTriangles topology;
  Points vertices( topology.points().dim() );

  Topology<Simplex<Lagrange>> topology_curved( vertices , topology , 3 );

  vertices.print();
  //builder.transfer( topology_curved );

}
UT_TEST_CASE_END( simplex_tests )

UT_TEST_SUITE_END( TopologySuite )
