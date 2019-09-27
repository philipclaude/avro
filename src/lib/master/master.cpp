#include "master/master.h"
#include "master/quadrature.h"

namespace ursa
{

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
