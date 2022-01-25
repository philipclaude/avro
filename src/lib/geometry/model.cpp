//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "geometry/model.h"

namespace avro
{

void
Model::get_entities( std::vector<Entity*>& entities ) const
{
  for (index_t k=0;k<body_.size();k++)
    body_[k]->get_entities(entities);
}

void
Model::get_tessellatable_entities( std::vector<Entity*>& entities ) const
{
  for (index_t k=0;k<body_.size();k++)
    body_[k]->get_tessellatable(entities);
}

} // avro
