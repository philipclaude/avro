//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_NUMERICS_FUNCTIONS_H_
#define avro_NUMERICS_FUNCTIONS_H_

#include "types.h"

namespace avro
{

namespace numerics
{

index_t factorial( const index_t i );
index_t binomial( const index_t n , const index_t k );
index_t nchoosek( index_t n , index_t k );
real_t  eps( const real_t& x );

} // numerics

} // avro

#endif
