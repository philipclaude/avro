//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_GEOMETRY_EGADS_MODEL_H_
#define avro_LIB_GEOMETRY_EGADS_MODEL_H_

#include "geometry/model.h"

#include "geometry/egads/egads_fwd.h"

#include <string>

namespace avro
{

namespace EGADS
{

class Body;
class Context;

class Model : public avro::Model
{
public:
  Model( Context* context );
  Model( Context* context , const std::string& filename , bool split=false );
  Model( const std::string& filename , bool split=false );
  Model( coord_t number );

  ~Model();

  void determine_number();

  const Context& context() const { return *context_; }

  Entity* find_entity( index_t id , int object_class ) const;

  void print() const;

private:
   Context* context_;
   ego object_;
   bool mine_;
};

} // EGADS

} // avro

#endif
