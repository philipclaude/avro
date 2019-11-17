#include "common/error.h"

#include "master/polytope.h"

#include "mesh/topology.h"

namespace luna
{

Polytope::Polytope( coord_t number , coord_t order , Data<int>& incidence ) :
  Master(number,order),
  simplex_(number,order),
  incidence_(incidence)
{
  luna_assert_msg( order==1 , "not supported..." );
}

Polytope::Polytope( Topology<Polytope>& topology , const coord_t order ) :
  Polytope(topology.number(),order,topology.master().incidence())
{}

} // luna
