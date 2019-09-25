#include "common/data.h"

#include "numerics/field.h"

namespace ursa
{

template<typename M,typename T>
Field<M,T>::Field( const Topology<M>& topology , int order ) :
  topology_(topology),
  master_(2,order)
{
  printf("constructing field with order %lu\n",master_.order());
}


template class Field< Master<LagrangeSimplex> , real >;

} // ursa
