#include "common/data.h"
#include "common/tree.hpp"

#include "mesh/topology.h"
#include "mesh/topology.hpp"

namespace ursa
{

/*\
 * Lagrange simplex topology
\*/
Topology<Simplex<Lagrange>>::Topology( Vertices& vertices , const Topology<Master_t>& linear , coord_t order ) :
  _Topology(vertices,linear.number(),order)
{
  convert(linear);
}

void
Topology<Simplex<Lagrange>>::convert( const Topology<Master_t>& linear )
{
  printf("converting to order %u...\n",master_.order());
}

/*\
 * Bezier simplex topology
\*/
Topology<Simplex<Bezier>>::Topology( Vertices& vertices , const Topology<Simplex<Lagrange>>& lagrange ) :
  _Topology(vertices,lagrange.number()),
  lagrange_(lagrange)
{
  convert();
}

void
Topology<Simplex<Bezier>>::convert()
{
  printf("convert lagrange basis of %u-simplices...\n",lagrange_.number());
  ursa_implement;
}

template class _Topology< Simplex<Lagrange> >;
template class Tree< Topology< Simplex<Lagrange> > >;

} // ursa
