//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "geometry/psc/node.h"

#include "numerics/geometry.h"

namespace avro
{

namespace PSC
{

void
Node::inverse( std::vector<real_t>& x , std::vector<real_t>& u ) const
{
  avro_assert( x.size() == dim_ );
  #if 0
  for (index_t j=0;j<x.size();j++)
    x[j] = (*this)(j);
  #else
  std::fill( u.begin() , u.end() , 1e20 );
  real_t d = numerics::distance( x_.data() ,  x.data() , dim_ );
  if (d < 1e-12)
  {
    u[0] = 0; 
  }
  else
    std::fill( x.begin() , x.end() , 1e20 );
    
  #endif
}



void
Node::evaluate( const std::vector<real_t>& , std::vector<real_t>& x ) const
{
  avro_assert( x.size() == dim_ );
  for (index_t j=0;j<x.size();j++)
    x[j] = (*this)(j);
}

} // PSC

} // avro
