#include "geometrics/egads.h"
#include "geometrics/plc.h"
#include "geometrics/primitive.h"
#include "geometrics/primitive.hpp"

#include "numerics/coordinate.h"

namespace ursa
{

namespace geometrics
{

template class Primitive<EGADS::Object>;
template class Primitive<PLC::Object>;

} // geometry

} // ursa
