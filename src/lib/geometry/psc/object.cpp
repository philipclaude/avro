#include "geometry/plc.h"

#include "numerics/coordinate.h"

namespace avro
{

namespace PLC
{

void
Object::inverse( numerics::Coordinate& x , numerics::Coordinate& u ) const
{
  x[0] = 0;
}

void
Object::evaluate( const numerics::Coordinate& u , numerics::Coordinate& x ) const
{
  x[0] = 0;
}

} // PLC

} // avro
