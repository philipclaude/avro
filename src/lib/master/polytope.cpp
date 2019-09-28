#include "common/error.h"

#include "master/polytope.h"

namespace ursa
{

Polytope::Polytope( coord_t number , coord_t order , Data<int>& incidence ) :
  Master(number,order),
  simplex_(number,order),
  incidence_(incidence)
{
  ursa_assert_msg( order==1 , "not supported..." );
}

} // ursa
