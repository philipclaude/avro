#include "master/master.h"
#include "master/quadrature.h"

namespace ursa
{

void
Simplex::precalculate()
{
}

void
Simplex::loadQuadrature( Quadrature& quadrature )
{
  quadrature.retrieve(xquad_,wquad_);
}

template class Master<LagrangeSimplex>;

} // ursa
