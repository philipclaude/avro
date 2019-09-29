#ifndef URSA_LIB_GEOMETRICS_PLC_H_
#define URSA_LIB_GEOMETRICS_PLC_H_

#include "common/types.h"

#include "numerics/matrix.h"

namespace ursa
{

namespace numerics
{
  class Coordinate;
}

namespace geometrics
{

namespace PLC
{

class Object
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

} // geometrics

} // ursa

#endif
