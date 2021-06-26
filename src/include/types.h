//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef AVRO_API_TYPES_H_
#define AVRO_API_TYPES_H_

#ifndef nil
#include <cstdlib>
#define nil NULL
#endif

namespace avro
{

typedef unsigned short coord_t;
typedef unsigned long  index_t;
typedef double         real_t;
typedef coord_t        coord_index_t;

enum Sign
{
  NEGATIVE = -1,
  ZERO = 0,
  POSITIVE = 1
};

} // avro

#endif
