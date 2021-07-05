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

#include "common/tools.h"

#include "geometry/egads/context.h"

#include "library/ckf.h"
#include "library/egads.h"
#include "library/tesseract.h"

#include "element/element.h"
#include "element/quadrature.h"

#include "mesh/topology.h"
#include "mesh/points.h"

using namespace avro;

UT_TEST_SUITE( mesh_curvilinear_suite )

UT_TEST_CASE( test1 )
{
  coord_t number = 4;
  coord_t dim = number;
  std::vector<index_t> dims(number,5);
  CKF_Triangulation linear_topology(dims);

  coord_t order = 3;
  Points points(dim);

  linear_topology.element().set_basis(BasisFunctionCategory_Lagrange);
  Topology<Simplex> curvilinear_topology(points,linear_topology,order);

  //points.print();
  //curvilinear_topology.Table<index_t>::print();
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( mesh_curvilinear_suite )
