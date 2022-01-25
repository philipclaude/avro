//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "unit_tester.hpp"

#include "adaptation/adapt.h"
#include "adaptation/geometry.h"
#include "adaptation/metric.h"

#include "common/error.h"

#include "geometry/egads/context.h"
#include "geometry/tessellation.h"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/egads.h"
#include "library/meshb.h"
#include "library/metric.h"
#include "library/tesseract.h"

#include "mesh/mesh.h"

using namespace avro;

UT_TEST_SUITE(insertions_geometry_suite)

UT_TEST_CASE(test1)
{
  // setup the topology
  coord_t number = 4;

  // parameters
  library::MetricField_Uniform analytic(number,0.1);

  // geometry
  std::vector<real_t> lens(number,1.0);
  std::vector<real_t> x0(number,0.5);
  library::Tesseract body(x0,lens);

  // structured grid
  std::vector<index_t> dims(number,3);
  CKF_Triangulation topology(dims);

  // tag the points onto the body
  topology.points().attach( body );

  topology.close();
  topology.neighbours().compute();
  topology.inverse().build();

  // create the metric field
  std::vector<symd<real_t>> fld;
  for (index_t k=0;k<topology.points().nb();k++)
    fld.push_back( analytic( topology.points()[k] ) );

  MetricAttachment attachment(topology.points(),fld);
  MetricField<Simplex> metric(topology,attachment);

  FieldInterpolation<Simplex,Metric> interpolation(&metric);
  metric.set_interpolation(&interpolation);
  topology.element().set_basis( BasisFunctionCategory_Lagrange );

  Smooth<Simplex> smoother(topology);
  UT_ASSERT( !smoother.element().parameter() );
  smoother.curved() = false;

  // loop through every edge, extract the cavity and ensure every edge end point is visible
  for (index_t iter=0;iter<10;iter++)
  for (index_t k=0;k<topology.points().nb();k++)
  {
    if (k < topology.points().nb_ghost()) continue;
    bool accept = smoother.apply( k , metric , -1 );
    UNUSED(accept);
  }

  //graphics::Viewer vis;
  //vis.add(topology);

  //vis.run(AVRO_FULL_UNIT_TEST);

}
UT_TEST_CASE_END(test1)


UT_TEST_SUITE_END(insertions_geometry_suite)
