#include "master/simplex.h"

namespace ursa
{

Simplex<Bezier>::Simplex( coord_t number , coord_t order ) :
  SimplexBase(number,order)
{
  // compute the transformation between lagrange and bezier bases

}


void
Simplex<Bezier>::evaluate( const real_t* x , real_t* y ) const
{
  ursa_implement;
}

void
Simplex<Bezier>::evaluate( index_t k , std::vector<real_t>& phi ) const
{
  // evaluate the set of basis functions at the k'th interpolation point

}

} // ursa
