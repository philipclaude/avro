#ifndef LUNA_LIB_GEOMETRY_PLC_H_
#define LUNA_LIB_GEOMETRY_PLC_H_

#include "common/types.h"

#include "geometry/entity.h"

#include "numerics/matrix.h"

namespace luna
{

namespace numerics
{
  class Coordinate;
}

namespace PLC
{

class Object : public Entity
{
public:
  Object( coord_t num , coord_t dim );

  void inverse( numerics::Coordinate& x , numerics::Coordinate& u ) const;
  void evaluate( const numerics::Coordinate& u , numerics::Coordinate& p ) const;

private:
  coord_t num_;
  coord_t dim_;

  numerics::MatrixD<real_t> basis_;
};

} // PLC

} // luna

#endif
