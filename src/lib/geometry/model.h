#ifndef luma_LIB_GEOMETRY_MODEL_H_
#define luma_LIB_GEOMETRY_MODEL_H_

#include "common/error.h"

#include "geometry/body.h"

namespace luma
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

protected:
  coord_t number_;

  std::vector<std::shared_ptr<Body>> body_;

};

} // luma

#endif
