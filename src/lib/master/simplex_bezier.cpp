#include "master/simplex.h"

namespace ursa
{

Simplex<Bezier>::Simplex( coord_t number , coord_t order ) :
  SimplexBase(number,order)
{
}


void
Simplex<Bezier>::evaluate( const real_t* x , real_t* y ) const
{
  ursa_implement;
}

} // ursa
