#include "geometry/egads.h"

#include "numerics/coordinate.h"

#include <egads.h>

namespace luna
{

namespace EGADS
{

Context::Context() :
  mine_(true)
{
  context_ = nil; // TODO create context
}

Context::Context( ego* context ) :
  context_(context),
  mine_(false)
{}

Context::~Context()
{
  if (mine_) delete context_;
}

ego*
Context::get()
{
  return context_;
}

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

Body::Body( const Context& context , ego* obj ) :
  luna::Body(0),
  context_(context),
  ego_(obj)
{
  // ask egads for the topological number
  number_ = 0;
}

void
Body::build()
{
  // call egads and then build the hierarchy
}

} // EGADS

} // luna
