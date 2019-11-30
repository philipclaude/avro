#include "common/tree.hpp"

#include "mesh/builder.h"
#include "mesh/topology.h"
#include "mesh/topology.hpp"

namespace luna
{

template<>
Topology<Simplex>::Topology( Points& vertices , coord_t number , coord_t order ) :
  TopologyBase(vertices,number,TableLayout_Rectangular),
  master_( number , order ),
  neighbours_(*this),
  inverse_(*this)
{}

template<>
Topology<Simplex>::Topology( Points& points , const Topology<Simplex>& linear , coord_t order ) :
 Topology(points,linear.number(),order)
{
  printf("converting to order %u...\n",master_.order());
  Builder<Simplex> builder(linear,master_.order(),BasisFunctionCategory_Lagrange);
  builder.transfer(*this);
}

template class Topology<Simplex>;
template class Tree<Topology<Simplex>>;

} // luna
