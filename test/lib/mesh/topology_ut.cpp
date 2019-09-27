#include "unit_tester.hpp"

#include "master/master.h"
#include "master/quadrature.h"

#include "mesh/topology.h"
#include "mesh/vertices.h"

using namespace ursa;

UT_TEST_SUITE( TopologySuite )

UT_TEST_CASE( simplex_tests )
{
  Vertices vertices(3);

  Topology< Simplex<Lagrange> > topology(vertices);

  topology.do_something();

}
UT_TEST_CASE_END( simplex_tests )

UT_TEST_SUITE_END( TopologySuite )
