#include "mesh/builder.h"

namespace ursa
{

template<typename Shape_t,typename Master_t>
template<typename T>
void
Builder<Shape_t,Master_t>::transfer( Field<Shape_t,Master_t,T>& F ) const
{
  ursa_implement;
}

template<typename Shape_t,typename Master_t>
template<typename MasterFrom_t,typename T>
void
Builder<Shape_t,Master_t>::transfer( const Field<Shape_t,MasterFrom_t,T>& fx , Field<Shape_t,Master_t,T>& fy ) const
{
  ursa_implement;
}

} // ursa
