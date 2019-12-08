#include "geometry/egads/object.h"

#include "numerics/coordinate.h"

#include <egads.h>

namespace luna
{

namespace EGADS
{

Object::Object( const Context& context , ego* obj ) :
  Entity(0),
  context_(context),
  ego_(obj)
{
  // get egads information to determine number...
}

ego*
Object::object()
{
  return ego_;
}

const ego*
Object::object() const
{
  return ego_;
}

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

} // EGADS

} // luna
