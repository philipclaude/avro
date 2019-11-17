#include "mesh/builder.h"

namespace luna
{

template<typename Shape_t,typename Master_t>
template<typename T>
void
Builder<Shape_t,Master_t>::transfer( Field<Shape_t,Master_t,T>& F ) const
{
  luna_implement;
}

template<typename Shape_t,typename Master_t>
template<typename MasterFrom_t,typename T>
void
Builder<Shape_t,Master_t>::transfer( const Field<Shape_t,MasterFrom_t,T>& fx , Field<Shape_t,Master_t,T>& fy ) const
{
  luna_implement;
}

} // luna
