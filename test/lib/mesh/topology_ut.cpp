#include "unit_tester.hpp"

#include "common/tools.h"

#include "master/master.h"
#include "master/quadrature.h"

#include "mesh/topology.h"
#include "mesh/vertices.h"

using namespace ursa;

UT_TEST_SUITE( TopologySuite )

UT_TEST_CASE( simplex_tests )
{
  Vertices vertices(3);
  coord_t number = 3;

  Topology< Simplex<Lagrange> > topology(vertices,number);

  topology.do_something();

  std::shared_ptr< Topology<Simplex<Lagrange>> > leaf = std::make_shared<Topology<Simplex<Lagrange>>>(vertices,number);
  topology.add_child(leaf);

  Topology<Simplex<Lagrange>>& c = topology.topology(0);
  UNUSED(c);

}
UT_TEST_CASE_END( simplex_tests )

UT_TEST_SUITE_END( TopologySuite )
