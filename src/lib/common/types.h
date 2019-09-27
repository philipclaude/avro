// ursa: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef URSA_COMMON_TYPES_H_
#define URSA_COMMON_TYPES_H_

namespace ursa
{

typedef unsigned short coord_t;
typedef unsigned long  index_t;
typedef double        real_t;

enum Sign {

    NEGATIVE = -1,

    ZERO = 0,

    POSITIVE = 1
};

typedef coord_t coord_index_t;

}

#ifndef nil
#include <cstdlib>
#define nil NULL
#endif

#endif
