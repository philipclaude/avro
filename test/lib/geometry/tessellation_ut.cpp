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

#include "geometry/entity.h"
#include "geometry/tessellation.h"

#include "geometry/egads/context.h"
#include "geometry/egads/model.h"

#include "graphics/application.h"

using namespace avro;

UT_TEST_SUITE( geometry_tessellation_suite )

UT_TEST_CASE(test1)
{
  EGADS::Context context;
  //EGADS::Model model(&context,BASE_TEST_DIR+"/geometry/tire.egads" );
  EGADS::Model model(&context,BASE_TEST_DIR+"/geometry/cube-cylinder.egads");
  //EGADS::Model model(&context,"/Users/pcaplan/Codes/EngSketchPad/data/basic/import_2.egads" );
  //EGADS::Model model(&context,"/Users/pcaplan/Codes/mach-II/library/geometry/turtle.step");

  TessellationParameters params;
  params.standard();

  params.min_size() = 0.1;
  params.min_angle() = 20;

  ModelTessellation tess(model,params);

  //tess.points().print(true);

  for (index_t k=0;k<tess.points().nb();k++)
  {
    if (tess.points().entity(k)->number()<=1)
      UT_ASSERT( tess.points().u(k,1) > 1e10 );
    else
      UT_ASSERT( tess.points().u(k,1) < 1e10 );
  }

  graphics::Visualizer vis;

  vis.add_topology( tess.topology(0) );

  vis.run();
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END( geometry_tessellation_suite )
