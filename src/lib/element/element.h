//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_SHAPE_SHAPE_H_
#define avro_LIB_SHAPE_SHAPE_H_

#include "avro_types.h"

#include <vector>

namespace avro
{

class Shape {
public:
  coord_t number() const { return number_; }
  coord_t order() const { return order_; }

protected:
  Shape( coord_t number , coord_t order ) :
    number_(number),
    order_(order)
  {}

  coord_t number_;
  coord_t order_;
};

typedef struct
{
  std::vector<index_t> indices;
  coord_t dim;
  bool sorted = true;
} ElementIndices;

bool operator< ( const ElementIndices& f , const ElementIndices& g );
bool operator== ( const ElementIndices& f , const ElementIndices& g );

} // avro

#endif
