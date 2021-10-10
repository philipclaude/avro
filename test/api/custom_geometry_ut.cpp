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
#include "numerics/geometry.h"

#include "avro.h"
#include "avro_config.h"

#include <egads.h>

#include <cmath>

using namespace avro;

UT_TEST_SUITE(api_custom_geometry_test_suite)

#if AVRO_NO_ESP == 0

void
add_children( EGADSGeneralGeometry& egg , Entity* entity ) {

  ego e = static_cast<EGADS::Object*>(entity)->object();
  for (index_t k = 0; k < entity->nb_children(); k++) {
    ego ek = static_cast<EGADS::Object*>(&entity->child(k))->object();
    egg.add_child(e,ek);
    add_children(egg,&entity->child(k));
  }
}

UT_TEST_CASE(test1)
{

  const coord_t dim = 3;
  coord_t number = dim;
  EGADSGeneralGeometry egg(number-1);

  // this example glues two boxes together and tags the interior face between them as an "interior"

  // create the first box
  ego box1;
  real_t params1[6] = { 0.0 , 0.0 , 0.0 , 1.0 , 1.0 , 0.5 };
  EGADS_ENSURE_SUCCESS( EG_makeSolidBody( egg.context()->get() , BOX , params1 , &box1 ) );

  // create the second box
  ego box2;
  real_t params2[6] = { 0.0 , 0.0 , 0.5 , 1.0 , 1.0 , 0.5 };
  EGADS_ENSURE_SUCCESS( EG_makeSolidBody( egg.context()->get() , BOX , params2 , &box2 ) );

  // create the avro bodies (this is only needed for this demo, ideally you would have your own wrapper around ego's)
  avro::EGADS::Body b1( *egg.context() , box1 );
  avro::EGADS::Body b2( *egg.context() , box2 );
  b1.build_hierarchy();
  b2.build_hierarchy();

  // add each body ego to the egg
  egg.add_body( box1 );
  egg.add_body( box2 );

  // retrieve the entities of the first body
  std::vector<Entity*> entities1;
  b1.get_entities(entities1);
  printf("body1 has %lu entities\n",entities1.size());

  // retrieve the entities of the second body
  std::vector<Entity*> entities2;
  b2.get_entities(entities2);
  printf("body2 has %lu entities\n",entities2.size());

  // add the ego's for each body (even if some are duplicates and won't be used)
  std::vector<Entity*> entities;
  for (index_t k = 0; k < entities1.size(); k++) {
    ego ek = static_cast<EGADS::Object*>(entities1[k])->object();
    egg.add_object( box1 , ek );
    entities.push_back( entities1[k] );
  }

  for (index_t k = 0; k < entities2.size(); k++) {
    ego ek = static_cast<EGADS::Object*>(entities2[k])->object();
    egg.add_object( box2 , ek );
    entities.push_back( entities2[k] );
  }
  printf("--> added objects\n");

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

  // attach the mesh points to the bodies (again, this is only for the demo)
  // if you created your own meshes to start
  // (say from TetGen, then you should retain which ego each vertex is on)
  mesh1.points().attach( b1 );
  mesh2.points().attach( b2 );

  // map duplicate entities (which should be at z = 0.5)
  // there should be 1 face, 4 edges and 4 nodes which are duplicates
  // (see the topological number printed below)
  std::map<Entity*,Entity*> entity_map;
  for (index_t i = 0; i < entities1.size(); i++)
  for (index_t j = 0; j < entities2.size(); j++) {

    ego ei = static_cast<EGADS::Object*>(entities1[i])->object();
    ego ej = static_cast<EGADS::Object*>(entities2[j])->object();

    int same = EG_isSame(ei,ej);
    if (same == EGADS_SUCCESS)
      printf("number = %u, i = %lu, j = %lu are the same!\n",entities1[i]->number(),i,j);
    else
      continue;

    // we will only keep entities1[i]
    entity_map.insert( {entities2[j],entities1[i]} );

    // mark the entity as interior
    entities1[i]->set_interior(true);

    // add the unique entity as a child of the duplicate entity's parents
    std::shared_ptr<Entity> ep = b1.lookup(ei); // retrieve the entity that we want to keep
    avro_assert( ep != nullptr );
    for (index_t k = 0; k < entities2[j]->nb_parents(); k++) {
      entities2[j]->parents(k)->add_child(ep);
    }
  }
  printf("--> number of duplicate entities = %lu\n",entity_map.size());

  // set the interior entities for the egg
  for (index_t k = 0; k < entities1.size(); k++) {
    if (!entities1[k]->interior()) continue;
    ego e = static_cast<EGADS::Object*>(entities1[k])->object();
    egg.set_interior(e);
  }

  // add the adjacency relationships
  for (index_t k = 0; k < entities.size(); k++) {
    add_children(egg,entities[k]);
  }
  printf("--> adjacencies added\n");

  // create one final mesh
  Points points(dim);
  Topology<Simplex> topology(points,number);

  // copy the points
  mesh1.points().copy( points );
  std::map<index_t,index_t> point_map;
  for (index_t k = 0; k < mesh2.points().nb(); k++) {

    Entity* ek = mesh2.points().entity(k);

    // check if the point is not a duplicate
    if (ek == nullptr || entity_map.find(ek) == entity_map.end()) {
      index_t idx = points.nb();
      point_map.insert( {k,idx} );
      points.create( mesh2.points()[k] );
      points.set_entity( idx , ek );
      points.set_param( idx , mesh2.points().u(k) );
      continue;
    }

    // determine which point in mesh1 this is closest to
    int idx = -1;
    for (index_t i = 0; i < mesh1.points().nb(); i++) {

      if (mesh1.points().entity(i) == nullptr) continue;
      if (!mesh1.points().entity(i)->interior()) continue; // reduce the loop down to O(# duplicate points)

      real_t d = numerics::distance( mesh1.points()[i] , mesh2.points()[k] , dim );
      if (d < 1e-12) {
        idx = i;
        break;
      }
    }
    avro_assert( idx >= 0 );
    point_map.insert( {k,idx} );
  }
  printf("--> number of duplicate points = %lu\n",
                     mesh1.points().nb() + mesh2.points().nb() - points.nb());
  UT_ASSERT_EQUALS( mesh1.points().nb() + mesh2.points().nb() - points.nb() , N*N );

  // compute the parameter coordinates (which may have changed when we adjusted for duplicates)
  for (index_t k = 0; k < points.nb(); k++) {
    Entity* e = points.entity(k);
    if (e == nullptr) continue;
    std::vector<real_t> X( points[k] , points[k] + dim );
    std::vector<real_t> U( dim-1 );
    e->inverse( X , U );
    points.set_param( k , U.data() );
  }

  // check the number of points on geometry Nodes is 12 (8 + 8 - the 4 duplicates)
  index_t nb_nodes = 0;
  for (index_t k = 0; k < points.nb(); k++) {
    if (points.entity(k) == nullptr) continue;
    if (points.entity(k)->number() == 0) nb_nodes++;
  }
  UT_ASSERT_EQUALS( nb_nodes , 12 );

  // copy the simplices
  for (index_t k = 0; k < mesh1.nb(); k++)
    topology.add( mesh1(k) , mesh1.nv(k) );
  for (index_t k = 0; k < mesh2.nb(); k++) {
    std::vector<index_t> elem( mesh2(k) , mesh2(k) + mesh2.nv(k) );
    for (index_t i = 0; i < elem.size(); i++)
      elem[i] = point_map.at( elem[i] );
    topology.add( elem.data() , elem.size() );
  }

  // finalize the customized geometry and pass it to the context
  egg.finalize();
  avro::Context adapter(dim,dim,dim-1);
  adapter.define_geometry(egg);

  // retrieve the association between ego's and the context's integer labels
  std::map<ego,int> ids;
  adapter.get_ego_ids(ids);
  printf("--> nb geometry ids = %lu\n",ids.size());

  // determine the geometry (integer) label for every vertex
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
    ego e = static_cast<EGADS::Object*>(entity)->object();
    geometry[k] = ids.at(e);
  }

  // send the mesh and geometry information to the context
  adapter.load_mesh( coordinates , topology.data() );
  adapter.load_geometry( geometry , parameters );

  index_t nb_points = points.nb();
  index_t nb_rank = number*(number+1)/2;
  std::vector<real_t> metrics( nb_points * nb_rank , 0.0 );

  // constant diagonal metric requesting a size of h
  real_t h = 0.15;
  for (index_t k = 0; k < nb_points; k++) {
    index_t idx = k*nb_rank;
    for (index_t i = 0; i < number; i++) {
      metrics[idx] = 1./(h*h);
      idx += number - i;
    }
  }

  //adapter.parameters().set( "output redirect" , "avro-output.txt" );
  adapter.parameters().set( "curved" , false );

  // perform the adaptation
  adapter.adapt(metrics);
}
UT_TEST_CASE_END(test1)

#endif // AVRO_NO_ESP

UT_TEST_SUITE_END(api_custom_geometry_test_suite)
