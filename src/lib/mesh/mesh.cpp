//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
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
