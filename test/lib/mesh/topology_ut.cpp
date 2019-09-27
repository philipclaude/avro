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

  std::shared_ptr< Topology<Simplex<Lagrange>> > leaf = std::make_shared<Topology<Simplex<Lagrange>>>(vertices);
  topology.addChild(leaf);

  Topology<Simplex<Lagrange>>& c = topology.child(0);

}
UT_TEST_CASE_END( simplex_tests )

UT_TEST_SUITE_END( TopologySuite )
