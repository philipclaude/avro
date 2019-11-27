#include "common/error.h"

#include "master/polytope.h"

#include "mesh/topology.h"

namespace luna
{

Polytope::Polytope( coord_t number , coord_t order ) :
  Master(number,order),
  simplex_(number,order)
{
  luna_assert_msg( order==1 , "not supported..." );
}

Polytope::Polytope( Topology<Polytope>& topology , const coord_t order ) :
  Master(topology.number(),order),
  simplex_(topology.number(),order)
{}


void
Polytope::get_edges( const index_t* v0 , index_t nv , std::vector<index_t>& ek ) const
{
  luna_implement;
}

} // luna
