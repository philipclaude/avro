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

#include "geometry/egads/context.h"
#include "geometry/egads/data.h"
#include "geometry/egads/object.h"

#include "graphics/application.h"

#include "library/ckf.h"

#include "avro.h"

#include <egads.h>

#include <cmath>

using namespace avro;

UT_TEST_SUITE(api_custom_geometry_test_suite)

UT_TEST_CASE(test1)
{
//#ifndef AVRO_NO_ESP

  const coord_t dim = 3;
  EGADSGeneralGeometry egg(dim-1);

  // this example glues two boxes together and tags the interior face between them as an "interior"
  EGADS::Context context;

  ego box1;
  real_t params1[6] = { 0.0 , 0.0 , 0.0 , 1.0 , 1.0 , 0.5 };
  EGADS_ENSURE_SUCCESS( EG_makeSolidBody( *context.get() , BOX , params1 , &box1 ) );

  ego box2;
  real_t params2[6] = { 0.0 , 0.0 , 0.5 , 1.0 , 1.0 , 0.5 };
  EGADS_ENSURE_SUCCESS( EG_makeSolidBody( *context.get() , BOX , params2 , &box2 ) );

  avro::EGADS::Body b1( context , box1 );
  avro::EGADS::Body b2( context , box2 );
  b1.build_hierarchy();
  b2.build_hierarchy();

  egg.add_body( box1 );
  egg.add_body( box2 );

  std::vector<Entity*> entities1;
  b1.get_entities(entities1);
  printf("body1 has %lu entities\n",entities1.size());

  std::vector<Entity*> entities2;
  b2.get_entities(entities2);
  printf("body2 has %lu entities\n",entities2.size());

  // add the ego's for each body (even if some are duplicates and won't be used)
  for (index_t k = 0; k < entities1.size(); k++) {
    ego* ek = static_cast<EGADS::Object*>(entities1[k])->object();
    egg.add_object( box1 , *ek );
  }

  for (index_t k = 0; k < entities2.size(); k++) {
    ego* ek = static_cast<EGADS::Object*>(entities2[k])->object();
    egg.add_object( box2 , *ek );
  }

  // generate a mesh in each box
  index_t N = 3;
  CKF_Triangulation mesh1( {N,N,N} );
  CKF_Triangulation mesh2( {N,N,N} );

  // squish the z-coordinates
  avro_assert( mesh1.points().nb() == mesh2.points().nb() );
  for (index_t k = 0; k < mesh1.points().nb(); k++) {
    mesh1.points()[k][2] *= 0.5;
    mesh2.points()[k][2] = mesh1.points()[k][2] + 0.5;
  }

  mesh1.points().attach( b1 );
  mesh2.points().attach( b2 );

  // look for the entity that is on z = 0.5
  std::map<Entity*,Entity*> entity_map;
  for (index_t i = 0; i < entities1.size(); i++)
  for (index_t j = 0; j < entities2.size(); j++) {

    ego ei = *static_cast<EGADS::Object*>(entities1[i])->object();
    ego ej = *static_cast<EGADS::Object*>(entities2[j])->object();

    int icode = EG_isSame(ei,ej);
    if (icode == EGADS_SUCCESS)
      printf("number = %u, i = %lu, j = %lu, is_same = %d\n",entities1[i]->number(),i,j,icode);

    entity_map.insert( {entities1[i],entities1[j]} );
    entities1[i]->set_interior(true);

    // add the unique entity as a child of the duplicate entity's parents
    std::shared_ptr<Entity> ep = b1.lookup(ei);
    for (index_t k = 0; k < entities2[j]->nb_parents(); k++) {
      entities2[j]->parents(k)->add_child(ep);
    }
  }

  for (index_t k = 0; k < mesh2.points().nb(); k++) {
    Entity* ek = mesh2.points().entity(k);
    if (ek == nullptr) continue;
    if (entity_map.find(ek) == entity_map.end()) continue;
    mesh2.points().set_entity(k,entity_map.at(ek));
  }

  // create one final mesh
  Points points(dim);
  Topology<Simplex> topology(points,dim);

  // copy the points
  mesh1.points().copy( points );
  for (index_t k = 0; k < mesh2.points().nb(); k++) {
    index_t idx = points.nb();
    points.create( mesh2.points()[k] );
    points.set_entity( idx , mesh2.points().entity(k) );
    points.set_param( idx , mesh2.points().u(k) );
  }

  for (index_t k = 0; k < mesh2.nb(); k++)
  for (index_t j = 0; j < mesh2.nv(k); j++)
    mesh2(k,j) += mesh1.points().nb();

  // copy the simplices
  for (index_t k = 0; k < mesh1.nb(); k++)
    topology.add( mesh1(k) , mesh1.nv(k) );
  for (index_t k = 0; k < mesh2.nb(); k++)
    topology.add( mesh2(k) , mesh2.nv(k) );

  // compute the parameter coordinates (which may have changed when we adjusted for duplicates)
  for (index_t k = 0; k < points.nb(); k++) {
    Entity* e = points.entity(k);
    if (e == nullptr) continue;
    std::vector<real_t> X( points[k] , points[k] + dim );
    std::vector<real_t> U( dim-1 );
    e->inverse( X , U );
    points.set_param( k , U.data() );
  }

  #if 1
  egg.finalize();
  avro::Context adapter(dim,dim,dim-1);
  adapter.define_geometry(egg);

  // retrieve the association between ego's and the context's integer labels
  std::map<int,int> ids;
  adapter.get_geometry_ids(ids);

  printf("nb geometry ids = %lu\n",ids.size());

  // define the mesh into the avro context
  std::vector<real_t> coordinates( dim*points.nb() );
  std::vector<int> geometry( points.nb() , -1 );
  std::vector<real_t> parameters( points.nb() * (dim-1) , unset_value );
  for (index_t k = 0; k < points.nb(); k++) {
    for (coord_t d = 0; d < dim; d++)
      coordinates[k*dim+d] = points[k][d];
    Entity* entity = points.entity(k);
    if (entity == nullptr) continue;
    parameters[2*k]   = points.u(k)[0];
    parameters[2*k+1] = points.u(k)[1];
    try {
      geometry[k] = ids.at(entity->identifier());
    }
    catch (...) {
      avro_assert_not_reached;
    }
  }
  adapter.load_mesh( coordinates , topology.data() );
  adapter.load_geometry( geometry , parameters );

  index_t nb_points = points.nb();
  coord_t number = dim;
  index_t nb_rank = number*(number+1)/2;
  std::vector<real_t> metrics( nb_points * nb_rank , 0.0 );

  // constant diagonal metric requesting a size of h
  real_t h = 0.25;
  for (index_t k = 0; k < nb_points; k++) {
    index_t idx = k*nb_rank;
    for (index_t i = 0; i < number; i++) {
      metrics[idx] = 1./(h*h);
      idx += number - i;
    }
  }

  points.print(true);

  adapter.parameters().set( "curved" , false );

  // perform the adaptation
  //adapter.adapt(metrics);

  #endif



  graphics::Viewer viewer;
  viewer.add( topology );
  viewer.run(AVRO_FULL_UNIT_TEST);

//#endif // AVRO_NO_ESP
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(api_custom_geometry_test_suite)
