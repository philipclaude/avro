#include "geometry/psc/node.h"

namespace avro
{

namespace PSC
{

void
Node::inverse( std::vector<real_t>& x , std::vector<real_t>& ) const
{
  avro_assert( x.size() == dim_ );
  for (index_t j=0;j<x.size();j++)
    x[j] = (*this)(j);
}



void
Node::evaluate( const std::vector<real_t>& , std::vector<real_t>& x ) const
{
  avro_assert( x.size() == dim_ );
  for (index_t j=0;j<x.size();j++)
    x[j] = (*this)(j);
}

} // PSC

} // avro
