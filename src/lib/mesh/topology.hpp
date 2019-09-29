#include "mesh/topology.h"
#include "mesh/vertices.h"

namespace ursa
{

template<typename type>
TopologyBase<type>::TopologyBase( Vertices& vertices , coord_t number , coord_t order ) :
  TopologyHolder(vertices,number),
  master_( number , order )
{
  printf("nb topologies = %lu\n",nb_children());
}

template<typename type>
TopologyBase<type>::TopologyBase( Vertices& vertices , coord_t number ) :
  TopologyBase(vertices,number,1)
{
  printf("nb topologies = %lu\n",nb_children());
}

} // ursa
