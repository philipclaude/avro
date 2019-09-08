#include "geometrics/egads.h"

#include "numerics/coordinate.h"

namespace ursa
{
namespace geometrics
{
namespace EGADS
{

Object::Object() :
  ego_(nullptr)
{

}

void
Object::project(numerics::Coordinate& x ) const
{
  x[0] = 0;
}

} // EGADS

} // geometrics

} // ursa
