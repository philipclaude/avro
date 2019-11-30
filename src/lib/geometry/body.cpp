#include "geometry/body.h"
#include "geometry/entity.h"

namespace luna
{

Body::Body( coord_t number ) :
  number_(number)
{}

void
Body::add( Entity_ptr entity )
{
  entity_.push_back(entity);
}

} // luna
