#include "common/error.h"

#include "master/polytope.h"

#include "mesh/topology.h"

namespace luma
{

Polytope::Polytope( coord_t number , coord_t order , const Table<int>& incidence ) :
  Master(number,order),
  simplex_(number,order),
  incidence_(incidence)
{
  luma_assert_msg( order==1 , "not supported..." );
}

Polytope::Polytope( Topology<Polytope>& topology , const coord_t order , const Table<int>& incidence ) :
  Master(topology.number(),order),
  simplex_(topology.number(),order),
  incidence_(incidence)
{}


void
Polytope::get_edges( const index_t* v0 , index_t nv , std::vector<index_t>& ek ) const
{
  luma_implement;
}

real_t
Polytope::volume( const Points& points , const index_t* v , index_t nv ) const
{
  luma_implement;
  return -1.0;
}


} // luma
