// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "unit_tester.hpp"

#include "geometry/entity.h"
#include "geometry/egads/context.h"

#include "library/ckf.h"
#include "library/egads.h"

#include "mesh/boundary.h"

using namespace avro;

UT_TEST_SUITE(boundary_test_suite)

#if 0
UT_TEST_CASE(test1)
{
  coord_t number = 2;

  std::vector<real> lens(number,1.);
  std::vector<index_t> dims(number,10);

  library::CubeMesh mesh( lens , dims );

  Context context;
  real_t xc[3] = {.5,.5,0.};
  library::EGADSSquare square( &context , xc , lens[0] , lens[1] );
  mesh.vertices().findGeometry(square);

  Topology<Simplex> topology(mesh.vertices(),number);
  mesh.retrieveElements(number,topology);
  Boundary<Simplex> boundary(topology);
  topology.neighbours().compute();
  boundary.extract();

  std::vector<Entity*> entities;
  square.listTessellatableEntities(entities);

  boundary.print();

  for (index_t k=0;k<entities.size();k++)
  {
    if (entities[k]->number()!=number-1) continue;
    index_t id = boundary.indexof( entities[k] );

    // calculate the volume
    real_t vol = boundary.volume( id );
    UT_ASSERT_NEAR( vol , 1. , 1e-12 );
  }
}
UT_TEST_CASE_END(test1)
#endif

UT_TEST_CASE(test2)
{
  EGADS::Context context;
  real_t xc[3] = {0.,0.,0.};

  bool full = false;

  EGADS::Cube box(&context,{1,1,1},xc);

  CKF_Triangulation topology({3,3,3});
  topology.points().attach(box);

  topology.points().print(true);

  Boundary<Simplex> boundary(topology);
  topology.neighbours().compute();
  if (full)
    boundary.extractall();
  else
    boundary.extract();

  if (full)
    UT_ASSERT_EQUALS( boundary.nb_children() , 26 );
  else
    UT_ASSERT_EQUALS( boundary.nb_children() , 6 );

  // retrieve all the entities so we can look up their associated topologies
  std::vector<Entity*> entities;
  box.get_tessellatable(entities);
  UT_ASSERT_EQUALS( entities.size() , 26 );

  for (index_t k=0;k<entities.size();k++)
  {
    if (entities[k]->number()!=2) continue;
    index_t id = boundary.indexof( entities[k] );

    // calculate the volume
    real_t vol = boundary.volume( id );
    UT_ASSERT_NEAR( vol , 1. , 1e-12 );
  }
}
UT_TEST_CASE_END(test2)

UT_TEST_CASE(test3)
{
  EGADS::Context context;
  real_t xc[3] = {0.,0.,0.};

  bool full = true;

  EGADS::Cube box(&context,{1,1,1},xc);

  CKF_Triangulation topology({3,3,3});
  topology.points().attach(box);

  topology.points().print(true);

  Boundary<Simplex> boundary(topology);
  topology.neighbours().compute();
  if (full)
    boundary.extractall();
  else
    boundary.extract();

  if (full)
    UT_ASSERT_EQUALS( boundary.nb_children() , 26 );
  else
    UT_ASSERT_EQUALS( boundary.nb_children() , 6 );

  // retrieve all the entities so we can look up their associated topologies
  std::vector<Entity*> entities;
  box.get_tessellatable(entities);
  UT_ASSERT_EQUALS( entities.size() , 26 );

  for (index_t k=0;k<entities.size();k++)
  {
    if (entities[k]->number()!=2) continue;
    index_t id = boundary.indexof( entities[k] );

    // calculate the volume
    real_t vol = boundary.volume( id );
    //UT_ASSERT_NEAR( vol , 1. , 1e-12 );
  }
}
UT_TEST_CASE_END(test3)



UT_TEST_SUITE_END(boundary_test_suite)
