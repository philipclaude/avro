#include "mesh/topology.h"
#include "mesh/points.h"

namespace avro
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
  set_layout( topology1.layout() );
  Table<index_t>::copy(topology1);
  closed_   = topology1.closed();
}

void
TopologyBase::set_points( Points& points )
{
  points_ = points;
  if (type_name()=="simplex")
    static_cast<Topology<Simplex>*>(this)->set_points(points);
  else
    avro_implement;
}

void
TopologyBase::offset_by( const index_t offset )
{
  // offset the indices
  for (index_t k=0;k<data_.size();k++)
    data_[k] += offset;
}

} // avro
