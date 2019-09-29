#include "common/tree.hpp"

#include "geometrics/egads.h"
#include "geometrics/plc.h"
#include "geometrics/primitive.h"

#include "numerics/coordinate.h"

namespace ursa
{

namespace geometrics
{

Primitive::Primitive( coord_t number ) :
  number_(number),
  name_("unnamed")
{}

Primitive::Primitive( coord_t number , const std::string& name ) :
  number_(number),
  name_(name)
{}

} // geometry

template class Tree<geometrics::Primitive>;

} // ursa
