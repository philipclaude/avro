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

#include "library/samples.h"

#include "mesh/builder.h"

using namespace avro;

UT_TEST_SUITE( builder_suite )

UT_TEST_CASE( simplex_tests )
{
  library::TwoTriangles topology;
  Points vertices( topology.points().dim() );

  topology.element().set_basis( BasisFunctionCategory_Lagrange );

  Topology<Simplex> topology_curved( vertices , topology , 2 );

  topology_curved.template Table<index_t>::print();

  vertices.print();
}
UT_TEST_CASE_END( simplex_tests )

UT_TEST_SUITE_END( builder_suite )
