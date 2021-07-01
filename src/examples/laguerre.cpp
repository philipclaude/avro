#include <avro.h>

#include <math.h>
#include <cstdlib>
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

  index_t nb_points = 100;

  std::vector<real_t> x( nb_points*dim );
  std::vector<real_t> w( nb_points , 0.0 );

  for (index_t k = 0; k < nb_points; k++)
  for (coord_t d = 0; d < dim; d++)
    x[dim*k+d] = real_t(rand())/real_t(RAND_MAX);

  context.compute_laguerre( x , w , 0 );

  std::vector<real_t> vertices;
  std::vector<index_t> polytopes;
  std::vector<index_t> nv_per_elem;
  context.retrieve_polytopes( vertices , polytopes , nv_per_elem );

  index_t nb_vertices = vertices.size()/dim;
  for (index_t k = 0; k < nb_vertices; k++) {
    printf("voronoi vertex[%lu] = ( ",k);
    for (coord_t d = 0; d < dim; d++)
      printf("%g ",vertices[k*dim+d]);
    printf(")\n");
  }

  index_t nb_polytopes = nv_per_elem.size();
  index_t i = 0;
  for (index_t k = 0; k < nb_polytopes; k++) {
    printf("polytope[%lu] = ( ",k);
    for (index_t j = 0; j < nv_per_elem[k]; j++)
      printf("%lu ",polytopes[i++]);
    printf(")\n");
  }

  assert( nb_polytopes == nb_points );

  return 0;
}
