#include "common/tools.h"

#include "mesh/topology.h"
#include "mesh/points.h"

#include <set>

namespace ursa
{

template<typename type>
TopologyBase<type>::TopologyBase( Points& vertices , coord_t number , coord_t order ) :
  TopologyHolder(vertices,number),
  master_( number , order )
{
  printf("nb topologies = %lu\n",nb_children());
}

template<typename type>
TopologyBase<type>::TopologyBase( Points& vertices , coord_t number ) :
  TopologyBase(vertices,number,1)
{
  printf("nb topologies = %lu\n",nb_children());
}

template<typename type>
void
TopologyBase<type>::getEdges( std::vector<index_t>& edges ) const
{
  std::vector<index_t> ek;

  std::set< std::pair<index_t,index_t> > table;

  for (index_t k=0;k<nb();k++)
  {
    if (ghost(k)) continue;

    const index_t* v0 = operator()(k);

    // get the edges of this cell
    master_.get_edges( v0 , nv(k) , ek );

    //printInline(ek,"ek");

    // add the edges
    for (index_t j=0;j<ek.size()/2;j++)
    {
      index_t p0 = ek[2*j];
      index_t p1 = ek[2*j+1];

      if (p0<points_.nb_ghost() || p1<points_.nb_ghost())
        continue;

      if (p0>p1) std::swap(p0,p1);
      std::pair<index_t,index_t> E = std::pair<index_t,index_t>(p0,p1);
      if (table.find(E)==table.end())
      {
        table.insert(E);
        edges.push_back(p0);
        edges.push_back(p1);
      }
    }
  }
}

} // ursa
