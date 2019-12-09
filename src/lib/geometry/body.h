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
protected:
  typedef std::shared_ptr<Entity> Entity_ptr;
public:
  virtual ~Body() {}
  
  coord_t number() const { return number_; }
  void add( Entity_ptr prim );

  index_t nb_entities() const { return entity_.size(); }

  void build_parents();

  virtual void print() const = 0;

protected:

  Body( coord_t number );

  std::vector<Entity_ptr> entity_;

  coord_t number_;
};

} // luna

#endif
