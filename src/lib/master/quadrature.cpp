#include "master/quadrature.h"

namespace ursa
{

Quadrature::Quadrature( coord_t dim , const int order) :
  dim_(dim),
  order_(order),
  defined_(false)
{}

void
ConicalProductQuadrature::define( const int order )
{
  printf("requesting order %d\n",order);

  // set nb_quad, x and w

  defined_ = true;
}

void
ConicalProductQuadrature::retrieve( std::vector<real>& x , std::vector<real>& w )
{
  if (!defined_) define(order_);


}



} // ursa
