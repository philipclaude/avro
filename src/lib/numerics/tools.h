// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef AVRO_NUMERICS_TOOLS_H_
#define AVRO_NUMERICS_TOOLS_H_

#include "common/types.h"

#include <vector>

namespace avro
{

namespace numerics
{

real_t eps( const real_t& x );
real_t sum( const std::vector<real_t>& x );
real_t exactsum( const std::vector<real_t>& x );
real_t naivesum( const std::vector<real_t>& x );
real_t norm( const std::vector<real_t>& x );

} // numerics

} // avro

#endif
