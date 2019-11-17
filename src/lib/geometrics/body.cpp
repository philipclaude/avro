#include "geometrics/body.h"
#include "geometrics/primitive.h"

namespace luna
{

namespace geometrics
{

Body::Body( coord_t number ) :
  number_(number)
{}

void
Body::add( Primitive_ptr prim )
{
  primitive_.push_back(prim);
}

} // geometrics

} // luna
