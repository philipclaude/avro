#ifndef LUNA_LIB_GEOMETRY_BODY_H_
#define LUNA_LIB_GEOMETRY_BODY_H_

#include "common/types.h"

#include <memory>
#include <vector>

namespace luna
{

class Entity;

class Body
{
  typedef std::shared_ptr<Entity> Entity_ptr;
public:
  coord_t number() const { return number_; }
  void add( Entity_ptr prim );

protected:

  Body( coord_t number );

  std::vector<Entity_ptr> entity_;

  coord_t number_;
};

} // luna

#endif
