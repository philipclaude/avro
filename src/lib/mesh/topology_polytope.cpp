#include "common/tree.hpp"

#include "mesh/builder.h"
#include "mesh/topology.h"
#include "mesh/topology.hpp"

namespace luna
{

template<>
Topology<Polytope>::Topology( Points& points , coord_t number , coord_t order ) :
  TopologyBase(points,number,TableLayout_Jagged),
  master_( number , order , points.incidence() ),
  neighbours_(*this),
  inverse_(*this)
{}

template class Topology<Polytope>;
template class Tree<Topology<Polytope>>;

} // luna
