#include "master/master.h"
#include "master/quadrature.h"

namespace ursa
{

template<typename Basis>
void
Simplex<Basis>::precalculate()
{
  ursa_implement;
}

template<typename Basis>
void
Simplex<Basis>::loadQuadrature( Quadrature& quadrature )
{
  quadrature.retrieve(xquad_,wquad_);
}

template class Master< Simplex<LagrangeSimplex> >;

} // ursa
