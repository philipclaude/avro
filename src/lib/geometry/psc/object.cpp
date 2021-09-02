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

#include "geometry/psc/object.h"

namespace avro
{

namespace PSC
{

Object::Object( coord_t number , coord_t dim ) :
  Entity(number),
  dim_(dim)
{
  tessellatable_ = true;
  interior_ = false;
}

void
Object::build_hierarchy()
{
  // not necesssary but declared virtual...maybe should just
  // use in egads geometry
  avro_assert_not_reached;
}

void
Object::print( bool with_children ) const
{
  printf("PSC %s entity with number %hu at %p\n",name_.c_str(),number_,(void*)(this));
  if (!with_children) return;
  for (index_t k=0;k<nb_children();k++)
    child(k).print(with_children);
}

} // PSC

} // avro
