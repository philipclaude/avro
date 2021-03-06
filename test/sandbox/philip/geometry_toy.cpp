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

#include "graphics/application.h"

#include "geometry/egads/context.h"
#include "geometry/entity.h"
#include "geometry/egads/model.h"
#include "geometry/tessellation.h"

#include "library/ckf.h"
#include "library/egads.h"

using namespace avro;

UT_TEST_SUITE( demo_toy )

UT_TEST_CASE( test1 )
{
  std::vector<index_t> dims(2,10);
  dims[1] = 5;

  CKF_Triangulation topology(dims);

  EGADS::Context context;
  std::vector<real_t> lens(2,1.);
  EGADS::Cube geometry(&context,lens);

  topology.points().attach(geometry);

  real_t volume = topology.volume();

  UT_ASSERT_NEAR( volume , 1.0 , 1e-12 );

  topology.points().print(true);

  index_t p = 29;
  Entity* entity = topology.points().entity(p);
  avro_assert( entity!=nullptr );

  Entity& entity_object = *entity;

  entity->print();
  //entity_object.print();

  EGADS::Model model(&context,"library/geometry/turtle.step");

  TessellationParameters params;
  params.set("min size" , 0.1 );
  params.set("min angle", 20.0);

  ModelTessellation tess(model,params);

  graphics::Viewer vis(true);

  vis.add( topology );
  vis.add( tess.topology(0) );

  vis.run();

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( demo_toy )
