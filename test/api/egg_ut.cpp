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
#include "graphics/math.h"

#include "library/ckf.h"
#include "numerics/geometry.h"
#include "numerics/linear_algebra.h"

#include "avro.h"
#include "avro_config.h"

#include <egads.h>

#include <cmath>

using namespace avro;

UT_TEST_SUITE(api_custom_geometry_test_suite)

#if AVRO_NO_ESP == 0

void
add_children( EGADSGeneralGeometry& egg , Entity* entity , std::map<ego,ego>& duplicates ) {

  ego e = static_cast<EGADS::Object*>(entity)->object();

  // do not add the duplicate ego's into the geometry hierarchy
  if (duplicates.find(e) != duplicates.end()) return;

  // add the children of this entity
  for (index_t k = 0; k < entity->nb_children(); k++) {
    ego ek = static_cast<EGADS::Object*>(&entity->child(k))->object();

    // again, do not add the duplicate ego's into the geometry hierarchy
    if (duplicates.find(ek) != duplicates.end()) continue;

    // add the child and recursively add its children as well
    egg.add_child(e,ek);
    add_children(egg,&entity->child(k),duplicates);
  }
}

UT_TEST_CASE(test1)
{
  // this example glues two boxes together and tags the interior Face between them as an "interior"
  // (note: interior Edges and Nodes on this Face are also tagged as interior)
  // please follow the different steps outlined below to set up your own custom geometry
  // (see the STEP ... and END STEP)

  // STEP 1: initialize  the geometry and avro context
  const coord_t dim = 3;
  coord_t number = dim;
  coord_t udim = dim -1;
  EGADSGeneralGeometry egg(number-1);
  avro::Context adapter(dim,dim,udim);
  // END STEP 1

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

  // STEP 2: add the ego's to the geometry
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
  // END STEP 2

  // generate a mesh in each box
  index_t N = 5;
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
  std::map<ego,ego> ego_map;
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
    ego_map.insert( {ej,ei} );

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

  // STEP 3: set the interior entities for the egg
  for (index_t k = 0; k < entities1.size(); k++) {
    if (!entities1[k]->interior()) continue;
    ego e = static_cast<EGADS::Object*>(entities1[k])->object();
    egg.set_interior(e);
  }
  // END STEP 3

  // STEP 4: add the ego parent-child relationships
  // add the adjacency relationships
  for (index_t k = 0; k < entities.size(); k++) {
    add_children(egg,entities[k],ego_map);
  }
  printf("--> adjacencies added\n");
  // END STEP 4

  // STEP 5: finalize the customized geometry and pass it to the context
  egg.finalize();
  adapter.define_geometry(egg);
  // END STEP 5

  // STEP 6: generate your mesh from the geometry
  // (in general you might use TetGen and retain the ego's each vertex is on)
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
  // END STEP 6

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

  // STEP 7: retrieve the association between ego's and the context's integer labels
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
  // END STEP 7, then proceed to use the context as usual

  // send the mesh and geometry information to the context
  adapter.load_mesh( coordinates , topology.data() );
  adapter.load_geometry( geometry , parameters );

  index_t nb_points = points.nb();
  index_t nb_rank = number*(number+1)/2;
  std::vector<real_t> metrics( nb_points * nb_rank , 0.0 );

  // constant diagonal metric requesting a size of h
  real_t h = 0.125;
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

  // retrieve the mesh from the context
  std::vector<index_t> elements;
  adapter.retrieve_mesh( coordinates , elements );
  nb_points = coordinates.size() / dim;
  index_t nb_elements = elements.size() / (number+1);
  printf("adapted mesh has %lu points and %lu elements\n",nb_points,nb_elements);

  // retrieve the geometry metadata
  std::vector<real_t> param;
  std::vector<int> entity_idx;
  adapter.retrieve_geometry( entity_idx , param );
  assert( param.size() == udim*nb_points );
  assert( entity_idx.size() == nb_points );
  for (index_t k = 0; k < nb_points; k++) {
    printf("\tpoint [%lu] (%g,%g,%g) on entity %d: u = ( ",k,coordinates[3*k],coordinates[3*k+1],coordinates[3*k+2],entity_idx[k]);
    for (index_t j = 0; j < udim; j++)
      printf("%g ",param[k*udim+j]);
    printf(")\n");
  }

  // retrieve the boundary facets
  std::vector<std::vector<index_t>> facets;

  // either of the following two functions will work (the first uses the current mesh stored in the context, the second uses the geometry information of each vertex along with the mesh elements to infer the boundary)
  // note that vertex coordinates are not required since this is a purely topological operation, which is why the second method does not require vertex coordinates
  //context.retrieve_boundary( facets , geometry );
  adapter.retrieve_boundary( entity_idx , elements , facets , geometry , true ); // true because we have interior entities

  assert( facets.size() == geometry.size() );
  std::vector<const real_t*> x(number);
  real_t area = 0.0;
  for (index_t bc = 0; bc < facets.size(); bc++) {
    const std::vector<index_t>& facets_on_bc = facets[bc];
    index_t nb_facets = facets_on_bc.size() / number;
    printf("there are %lu facets on bc %lu with geometry id %d\n",nb_facets,bc,geometry[bc]);

    for (index_t k = 0; k < nb_facets; k++) {


      // calculate the area of the facet
      for (index_t j = 0; j < number; j++)
        x[j] = &coordinates[3*facets_on_bc[k*number+j]];

      area += numerics::volume_nd( x , dim );
    }
  }

  // we have 2 external 1x1 faces, 8 external 0.5x1 faces, 1 internal 1x1 face
  printf("area = %g\n",area);
  UT_ASSERT_NEAR( area , 7.0 , 1e-12 );

}
UT_TEST_CASE_END(test1)

#endif // AVRO_NO_ESP

UT_TEST_SUITE_END(api_custom_geometry_test_suite)
