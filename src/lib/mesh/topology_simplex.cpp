#include "common/data.h"
#include "common/tree.hpp"

#include "mesh/topology.h"
#include "mesh/topology.hpp"

namespace ursa
{

template<typename Basis>
Topology< Simplex<Basis> >::Topology( Vertices& vertices ) :
  _Topology< Simplex<Basis> >(vertices,vertices.dim())
{
  printf("i'm a topology for %lu-simplices\n",this->master_.number());
}

template class Topology< Simplex<Lagrange> >;


} // ursa
