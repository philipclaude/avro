#include "common/data.h"
#include "common/tree.hpp"

#include "mesh/builder.h"
#include "mesh/topology.h"
#include "mesh/topology.hpp"

namespace luna
{

Topology<Simplex>::Topology( Points& points , coord_t number ) :
  TopologyBase<Simplex>(points,number)
{}

Topology<Simplex>::Topology( Points& points , const Topology<Simplex>& linear , coord_t order ) :
 TopologyBase(points,linear.number(),order)
{
  convert(linear);
}

void
Topology<Simplex>::convert( const Topology<Simplex>& linear )
{
  printf("converting to order %u...\n",master_.order());
  Builder<Simplex> builder(linear,master_.order(),BasisFunctionCategory_Lagrange);
  builder.transfer(*this);
}

template class TopologyBase<Simplex>;
template class Tree<Topology<Simplex>>;

} // luna
