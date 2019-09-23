#include "master/master.h"
#include "master/quadrature.h"

namespace ursa
{

template<typename Basis>
void
Simplex<Basis>::precalculate()
{
}

template<typename Basis>
void
Simplex<Basis>::loadQuadrature( Quadrature& quadrature )
{
  quadrature.retrieve(xquad_,wquad_);
}

template class Simplex<LagrangeSimplex>;

} // ursa
