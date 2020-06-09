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

#include "geometry/egads/context.h"
#include "geometry/egads/model.h"

#include "graphics/application.h"

#include "library/meshb.h"

using namespace avro;

UT_TEST_SUITE( meshb_test_suite )

UT_TEST_CASE( test1 )
{
  EGADS::Context context;
  EGADS::Model model(&context,BASE_TEST_DIR+"/geometry/cube-cylinder.egads");

  library::meshb mesh( BASE_TEST_DIR+"/meshes/cube-cylinder.mesh" , &model );

  // the original file does not have VerticesOn[entity] information, so we need to attach
  mesh.points().attach(model);
  mesh.points().compute_params();

  mesh.write( mesh , "tmp/cc1.mesh" , true );
  library::meshb mesh_in( "tmp/cc1.mesh" , &model );
  UT_ASSERT_EQUALS( mesh_in.nb_topologies() , 8 ); // 7 faces + 1 volume

  mesh.write( mesh_in , "tmp/cc2.mesh" , true );
  library::meshb mesh_in2( "tmp/cc2.mesh" , &model );
  UT_ASSERT_EQUALS( mesh_in2.nb_topologies() , 8 ); // 7 faces + 1 volume

  mesh.write( mesh_in2 , "tmp/cc3.mesh" , true );

  graphics::Visualizer vis;

  for (index_t k=0;k<mesh_in.nb_topologies();k++)
    vis.add_topology(mesh_in.topology(k));

  vis.run();
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( meshb_test_suite )
