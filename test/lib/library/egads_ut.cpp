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
#include "geometry/tessellation.h"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/egads.h"

#include "mesh/boundary.h"

using namespace avro;

UT_TEST_SUITE( library_egads_suite )

UT_TEST_CASE( box_test )
{
  EGADS::Context context;
  real_t x0[3] = {0,0,0};
  EGADS::Cube cube(&context,{1,1,1},x0);
}
UT_TEST_CASE_END( box_test )

UT_TEST_CASE( square_test )
{
  EGADS::Context context;
  EGADS::Cube square(&context,{1,1});

  CKF_Triangulation ckf( {10,10} );
  ckf.points().attach( square );

  Boundary<Simplex> boundary(ckf);
  boundary.extract();

  ckf.points().print(true);
}
UT_TEST_CASE_END( square_test )

UT_TEST_CASE( solid_bodies )
{
  EGADS::Context context;

  // sphere
  real_t rs = 1.0;
  real_t x0[3] = {0,0,0};
  std::shared_ptr<Body> sphere = std::make_shared<EGADS::Sphere>(&context,x0,rs);

  // cone
  real_t apex[3] = {0,0,1};
  std::shared_ptr<Body> cone = std::make_shared<EGADS::Cone>(&context,apex,x0,rs);

  // torus
  

  EGADS::Model model(2);
  model.add_body( sphere );
  model.add_body( cone );

  TessellationParameters params;
  params.standard();

  ModelTessellation tess(model,params);

  graphics::Visualizer vis;
  for (index_t j=0;j<tess.nb_topologies();j++)
    vis.add_topology(tess.topology(j));

  vis.run();


}
UT_TEST_CASE_END( solid_bodies )

UT_TEST_SUITE_END( library_egads_suite )
