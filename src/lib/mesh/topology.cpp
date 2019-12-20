#include "mesh/topology.h"
#include "mesh/points.h"

namespace luna
{

Fields::Fields( const TopologyBase& topology ) :
  topology_(topology)
{}

void
TopologyBase::copy( const TopologyBase& topology1 )
{
  dummy_    = topology1.dummy();
  number_   = topology1.number();
  name_     = topology1.name();
  Table<index_t>::copy(topology1);
  closed_   = topology1.closed();
}

} // luna
