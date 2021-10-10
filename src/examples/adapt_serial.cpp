//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include <avro.h>

#include <assert.h>

int
main( int argc , char** argv )
{
  using avro::index_t;
  using avro::coord_t;
  using avro::real_t;

  // setup the context
  coord_t number = 3;
  coord_t dim = number;
  coord_t udim = dim-1;
  avro::Context context(number,dim,udim);

  // setup the geometry & initial mesh
  context.define_geometry("box"); //
  context.define_mesh("CKF-3-3-3");
  context.attach_geometry();

  // retrieve the mesh coordinates so we know how many metrics to define
  std::vector<real_t> coordinates;
  std::vector<index_t> elements; // unused
  context.retrieve_mesh( coordinates , elements );
  index_t nb_points = coordinates.size() / dim;

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

  // overwrite some of the default parameters
  context.parameters().set( "output redirect" , "avro-output.txt" );

  // perform the adaptation
  context.adapt(metrics);

  // retrieve the mesh from the context
  context.retrieve_mesh( coordinates , elements );
  nb_points = coordinates.size() / dim;
  index_t nb_elements = elements.size() / (number+1);
  printf("adapted mesh has %lu points and %lu elements\n",nb_points,nb_elements);

  // retrieve the geometry metadata
  std::vector<real_t> param;
  std::vector<int> entity_idx;
  context.retrieve_geometry( entity_idx , param );
  assert( param.size() == udim*nb_points );
  assert( entity_idx.size() == nb_points );
  for (index_t k = 0; k < nb_points; k++) {
    printf("\tpoint [%lu] on entity %d: u = ( ",k,entity_idx[k]);
    for (index_t j = 0; j < udim; j++)
      printf("%g ",param[k*udim+j]);
    printf(")\n");
  }

  // retrieve the boundary facets
  std::vector<std::vector<index_t>> facets;
  std::vector<int> geometry;

  // either of the following two functions will work (the first uses the current mesh stored in the context, the second uses the geometry information of each vertex along with the mesh elements to infer the boundary)
  // note that vertex coordinates are not required since this is a purely topological operation, which is why the second method does not require vertex coordinates
  //context.retrieve_boundary( facets , geometry );
  context.retrieve_boundary( entity_idx , elements , facets , geometry );

  assert( facets.size() == geometry.size() );
  for (index_t bc = 0; bc < facets.size(); bc++) {
    const std::vector<index_t>& facets_on_bc = facets[bc];
    index_t nb_facets = facets_on_bc.size() / number;
    printf("there are %lu facets on bc %lu with geometry id %d\n",nb_facets,bc,geometry[bc]);

    for (index_t k = 0; k < nb_facets; k++) {
      printf("\tfacet[%lu] = ( ",k);
      for (index_t j = 0; j < number; j++)
        printf("%lu ",facets_on_bc[k*number+j]);
      printf(")\n");
    }
  }

  return 0;
}
