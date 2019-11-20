#include "mesh/builder.h"

namespace luna
{

template<typename type>
template<typename T>
void
Builder<type>::transfer( Field<type,T>& F ) const
{
  luna_implement;
}

template<typename type>
template<typename T>
void
Builder<type>::transfer( const Field<type,T>& fx , Field<type,T>& fy ) const
{
  luna_implement;
}

} // luna
