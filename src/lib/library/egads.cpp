//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "geometry/egads/context.h"
#include "library/egads.h"

#include <array>

namespace avro
{

namespace EGADS
{

Cube::Cube( const Context* context , coord_t number ) :
  Body(*context,&object_)
{}

Cube::Cube( const Context* context , const std::vector<real_t>& lens , const real_t* x0 ) :
  Cube(context,lens.size())
{
  std::vector<real_t> x(3,0);
  if (x0!=nullptr)
    x.assign(x0,x0+3);
  if (lens.size()==2)
  {
    real_t data[6] = { x[0] , x[1] , x[2] , lens[0] , lens[1] , lens[2] };
    EGADS_ENSURE_SUCCESS( EG_makeSolidBody( *context_.get() , BOX , data , &object_ ) );
  }
  else if (lens.size()==3)
  {
    real_t data[6] = { x[0] , x[1] , x[2] , lens[0] , lens[1] , lens[2] };
    EGADS_ENSURE_SUCCESS( EG_makeSolidBody( *context_.get() , BOX , data , &object_ ) );
  }
  this->build_hierarchy();
}

} // EGADS

} // avro
