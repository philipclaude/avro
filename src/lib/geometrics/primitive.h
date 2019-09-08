#ifndef URSA_LIB_GEOMETRICS_PRIMITIVE_H_
#define URSA_LIB_GEOMETRICS_PRIMITIVE_H_

#include "common/tree.h"

namespace ursa
{

namespace numerics
{
  class Coordinate;
}

namespace geometrics
{

class PrimitiveBase
{
protected:
  virtual ~PrimitiveBase() {}
  virtual void project( numerics::Coordinate& x ) const = 0;
};

template<typename T>
class Primitive : public Tree<Primitive<T>>, public PrimitiveBase
{
public:
  Primitive( T& prim ) :
    prim_(prim)
  {}

  ~Primitive() {}

  void project( numerics::Coordinate& x ) const;

private:
  T& prim_;
};

} // geometrics

} // ursa

#endif
