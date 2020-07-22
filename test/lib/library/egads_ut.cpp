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

#ifndef AVRO_NO_ESP
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
  real_t dir[3] = {0,0,1};
  real_t r0 = 0.1;
  real_t r1 = 2.0;
  std::shared_ptr<Body> torus = std::make_shared<EGADS::Torus>(&context,x0,dir,r1,r0);

  // cylinder
  real_t x1[3] = {0,0,1};
  real_t rc = 0.5;
  std::shared_ptr<Body> cylinder = std::make_shared<EGADS::Cylinder>(&context,x0,x1,rc);

  EGADS::Model model(2);
  model.add_body( sphere );
  model.add_body( cone );
  model.add_body( torus );
  model.add_body( cylinder );

  TessellationParameters params;
  params.standard();

  ModelTessellation tess(model,params);

  graphics::Visualizer vis;
  for (index_t j=0;j<tess.nb_topologies();j++)
    vis.add_topology(tess.topology(j));

  //vis.run();

}
UT_TEST_CASE_END( solid_bodies )

UT_TEST_CASE( smiley_test )
{
  EGADS::Context context;
  real_t x0[3] = {0,0,0};
  std::shared_ptr<Body> smiley = std::make_shared<EGADS::Smiley>( &context , x0 , 1.0 , 0.5 , 0.1 , M_PI/8 , 0.2 , 0.45 , M_PI/5. );

  EGADS::Model model(2);
  model.add_body(smiley);

  TessellationParameters params;
  params.standard();
  params.min_size() = 0.1;

  ModelTessellation tess(model,params);

  graphics::Visualizer vis;
  for (index_t j=0;j<tess.nb_topologies();j++)
    vis.add_topology(tess.topology(j));
  vis.run();
}
UT_TEST_CASE_END( smiley_test )

#endif

UT_TEST_SUITE_END( library_egads_suite )
