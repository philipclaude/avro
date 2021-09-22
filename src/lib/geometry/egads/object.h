//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_GEOMETRY_EGADS_H_
#define avro_LIB_GEOMETRY_EGADS_H_

#include "geometry/egads/body.h"
#include "geometry/egads/data.h"
#include "geometry/entity.h"

struct egObject;
typedef egObject* ego;

namespace avro
{

namespace numerics
{
  class Coordinate;
}

namespace EGADS
{

class Context;

class Object : public Entity
{
public:
  Object( const Context& context , ego* object );
  Object( ego* object , EGADS::Body* body );
  Object( const Context& context );

  void inverse( std::vector<real_t>& x , std::vector<real_t>& u ) const;
  void inverse_guess( std::vector<real_t>& x , std::vector<real_t>& u ) const;
  void evaluate( const std::vector<real_t>& u , std::vector<real_t>& p ) const;

  void project( std::vector<real_t>& x , std::vector<real_t>& u ) const;

  void set_object( ego* object );
  void construct( ego* object );

  void delete_object();

  void build_hierarchy();

  ego* object();
  ego* object() const;

  ego egchild( index_t k ) const;

  int object_class() const { return data_.object_class; }
  int member_type() const { return data_.member_type; }


  void print(bool with_children=true) const;

protected:
  EGADS::Body* body_;
  const Context& context_;
  ego* object_;

  egoData data_;

  int body_index_;
};

} // EGADS

} // avro

#endif
