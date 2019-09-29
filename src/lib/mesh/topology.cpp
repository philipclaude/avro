#include "mesh/topology.h"
#include "mesh/vertices.h"

namespace ursa
{

void
TopologyHolder::copy( TopologyHolder& topology1 )
{
  dummy_    = topology1.dummy();
  number_   = topology1.number();
  name_     = topology1.name();
  elements_ = topology1.elements();
  first_    = topology1.first();
  last_     = topology1.last();
  vertices_ = topology1.vertices();
}

} // avro
