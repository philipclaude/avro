//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "common/tools.h"

#include <math.h>
#include <cstdlib>

namespace avro
{

  real_t
  random_within( const real_t lo , const real_t hi )
  {
  	real_t s = (real_t) rand() / (real_t) RAND_MAX;
    return lo +s*(hi-lo);
  }

  int
  random_within( const int lo , const int hi )
  {
  	real_t s = (real_t)rand()/(real_t) RAND_MAX;
  	return lo +floor(s)*(hi -lo);
  }


}
