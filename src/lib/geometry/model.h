#ifndef avro_LIB_GEOMETRY_MODEL_H_
#define avro_LIB_GEOMETRY_MODEL_H_

#include "common/error.h"

#include "geometry/body.h"

namespace avro
{

class Entity;

class Model
{
public:
  Model( coord_t number ) :
    number_(number)
  {}

  Body& body(index_t k) { return *body_[k].get(); }
  const Body& body(index_t k) const { return *body_[k].get(); }

  index_t nb_bodies() const { return body_.size(); }

  void get_entities( std::vector<Entity*>& entities ) const;

  void add_body( std::shared_ptr<Body>& body )
  { body_.push_back(body); }

protected:
  coord_t number_;

  std::vector<std::shared_ptr<Body>> body_;

};

} // avro

#endif
