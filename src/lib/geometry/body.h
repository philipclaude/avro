//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_GEOMETRY_BODY_H_
#define avro_LIB_GEOMETRY_BODY_H_

#include "common/types.h"

#include <memory>
#include <vector>

namespace avro
{

class Entity;
class BodyTessellation;
class TessellationParameters;

class Body
{
protected:
  typedef std::shared_ptr<Entity> Entity_ptr;

public:

  coord_t number() const { return number_; }
  void add( Entity_ptr prim );

  index_t nb_entities() const { return entity_.size(); }

  void build_parents();

  virtual void print() const = 0;

  void get_entities( std::vector<Entity*>& entities ) const;
  void get_tessellatable( std::vector<Entity*>& entities ) const;

  virtual void tessellate( BodyTessellation& body_tess ) const = 0;

protected:

  Body( coord_t number );
  virtual ~Body() {}

  std::vector<Entity_ptr> entity_;

  coord_t number_;
};

} // avro

#endif
