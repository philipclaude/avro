#include <avro.h>

#include "common/process.h" // to initialize MPI

#define VERBOSE 0

int
main( int argc , char** argv )
{
  using avro::index_t;
  using avro::coord_t;
  using avro::real_t;

  avro::ProcessMPI::initialize(); // initialize MPI how you would in your code
  index_t rank = avro::ProcessMPI::rank();

  // setup the context
  coord_t number = 3;
  coord_t dim = number;
  coord_t udim = dim-1;
  avro::Context context(number,dim,udim);

  // setup the geometry & initial mesh
  context.define_geometry("box");
  context.define_mesh("CKF-3-3-3");
  context.attach_geometry();

  // the mesh is not partitioned yet, so avro will partition it with parmetis
  // you may not need to do something like this if you partition the mesh yourself
  // upon completion, the local mesh will be stored in the context
  context.partition();

  // retrieve the mesh coordinates so we know how many metrics to define
  std::vector<real_t> coordinates;
  std::vector<index_t> elements; // unused
  context.retrieve_mesh( coordinates , elements );
  index_t nb_points = coordinates.size() / dim; // local number of points

  index_t nb_rank = number*(number+1)/2;
  std::vector<real_t> metrics( nb_points * nb_rank , 0.0 ); // local metrics

  // constant diagonal metric requesting a size of h
  real_t h = 0.25;
  for (index_t k = 0; k < nb_points; k++) {
    index_t idx = k*nb_rank;
    for (index_t i = 0; i < number; i++) {
      metrics[idx] = 1./(h*h);
      idx += number - i;
    }
  }

  // perform the adaptation
  context.parameters().set( "limit metric" , true );
  context.adapt_parallel(metrics);

  // retrieve the mesh from the context
  context.retrieve_mesh( coordinates , elements );
  nb_points = coordinates.size() / dim;
  index_t nb_elements = elements.size() / (number+1);
  printf("[rank %lu]: adapted mesh has %lu points and %lu elements\n",rank,nb_points,nb_elements);

  // retrieve the geometry metadata
  std::vector<real_t> param;
  std::vector<int> entity_idx;
  context.retrieve_geometry( entity_idx , param );
  assert( param.size() == udim*nb_points );
  assert( entity_idx.size() == nb_points );
  #if VERBOSE
  for (index_t k = 0; k < nb_points; k++) {
    printf("\tpoint [%lu] on entity %d: u = ( ",k,entity_idx[k]);
    for (index_t j = 0; j < udim; j++)
      printf("%g ",param[k*udim+j]);
    printf(")\n");
  }
  #endif

  // retrieve the boundary facets
  std::vector<std::vector<index_t>> facets;
  std::vector<int> geometry;
  context.retrieve_boundary_parallel( facets , geometry );
  assert( facets.size() == geometry.size() );
  for (index_t bc = 0; bc < facets.size(); bc++) {
    const std::vector<index_t>& facets_on_bc = facets[bc];
    index_t nb_facets = facets_on_bc.size() / number;
    printf("there are %lu facets on bc %lu with geometry id %d\n",nb_facets,bc,geometry[bc]);

    #if VERBOSE
    for (index_t k = 0; k < nb_facets; k++) {
      printf("\tfacet[%lu] = ( ",k);
      for (index_t j = 0; j < number; j++)
        printf("%lu ",facets_on_bc[k*number+j]);
      printf(")\n");
    }
    #endif
  }

  avro::ProcessMPI::terminate();
}
