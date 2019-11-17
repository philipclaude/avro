#ifndef LUNA_LIB_GEOMETRICS_BODY_H_
#define LUNA_LIB_GEOMETRICS_BODY_H_

#include "common/types.h"

#include <memory>
#include <vector>

namespace luna
{

namespace geometrics
{

class Primitive;

class Body
{
  typedef std::shared_ptr<Primitive> Primitive_ptr;
public:
  coord_t number() const { return number_; }
  void add( Primitive_ptr prim );

protected:

  Body( coord_t number );

  std::vector<Primitive_ptr> primitive_;

  coord_t number_;
};

} // geometrics

} // luna

#endif
