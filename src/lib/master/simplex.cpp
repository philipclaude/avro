#include "master/quadrature.h"
#include "master/simplex.h"

#include "mesh/topology.h"

namespace ursa
{

template<typename Basis>
SimplexBase<Basis>::SimplexBase( const Topology<Simplex<Basis>>& topology , const coord_t order ) :
  SimplexBase(topology.number(),order)
{}

template<typename Basis>
void
SimplexBase<Basis>::precalculate()
{
  //ursa_implement;
}

template<typename Basis>
void
SimplexBase<Basis>::loadQuadrature( Quadrature& quadrature )
{
  quadrature.retrieve(xquad_,wquad_);
}

template class SimplexBase<Lagrange>;
template class SimplexBase<Bezier>;


} // ursa
