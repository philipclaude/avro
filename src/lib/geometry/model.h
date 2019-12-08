#ifndef LUNA_LIB_GEOMETRY_MODEL_H_
#define LUNA_LIB_GEOMETRY_MODEL_H_

#include "common/error.h"

#include "geometry/body.h"

namespace luna
{

class Model
{
protected:
  typedef std::shared_ptr<Body> Body_ptr;

public:
  Model( coord_t number ) :
    number_(number)
  {}

  Body& body(index_t k) { return *body_[k].get(); }
  const Body& body(index_t k) const { return *body_[k].get(); }

  void determine_number() { luna_implement; }

protected:
  coord_t number_;
  std::vector<Body_ptr> body_;
};

} // luna

#endif
