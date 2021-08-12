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

#include "geometry/body.h"
#include "geometry/entity.h"

namespace avro
{

Body::Body( coord_t number ) :
  number_(number)
{}

void
Body::add( Entity_ptr entity )
{
  entity_.push_back(entity);
}

void
Body::build_parents()
{
  for (index_t k=0;k<nb_entities();k++) {
    entity_[k]->build_parents();
    entity_[k]->build_children();
  }
}

void
Body::get_entities( std::vector<Entity*>& entities ) const
{
  for (index_t k=0;k<entity_.size();k++)
  {
    entities.push_back(entity_[k].get());
    entity_[k]->get_children(entities);
  }
  uniquify(entities);
}

void
Body::get_tessellatable( std::vector<Entity*>& entities ) const
{
  std::vector<Entity*> entities0;
  get_entities(entities0);
  for (index_t k=0;k<entities0.size();k++)
  {
    if (entities0[k]->tessellatable())
      entities.push_back(entities0[k]);
  }
}

} // avro
