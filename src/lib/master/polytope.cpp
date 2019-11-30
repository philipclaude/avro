#include "common/error.h"

#include "master/polytope.h"

#include "mesh/topology.h"

namespace luna
{

Polytope::Polytope( coord_t number , coord_t order , const Table<int>& incidence ) :
  Master(number,order),
  simplex_(number,order),
  incidence_(incidence)
{
  luna_assert_msg( order==1 , "not supported..." );
}

Polytope::Polytope( Topology<Polytope>& topology , const coord_t order , const Table<int>& incidence ) :
  Master(topology.number(),order),
  simplex_(topology.number(),order),
  incidence_(incidence)
{}


void
Polytope::get_edges( const index_t* v0 , index_t nv , std::vector<index_t>& ek ) const
{
  luna_implement;
}

} // luna
