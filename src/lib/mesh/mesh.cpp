#include "mesh/mesh.h"

namespace avro
{

Mesh::Mesh( coord_t dim ) :
  points_(dim),
  number_(dim)
{}

Mesh::Mesh( coord_t number , coord_t dim ) :
  points_(dim),
  number_(number)
{}

} // avro
