#include "mesh/builder.h"

namespace luma
{

template<typename type>
template<typename T>
void
Builder<type>::transfer( Field<type,T>& F ) const
{
  // transfer the associativity
  for (index_t k=0;k<nb();k++)
    F.add( (*this)(k) , nv(k) );

  // transfer the dof
  const std::vector<index_t>& idx = Table<index_t>::data();
  index_t nb_dof = * std::max_element( idx.begin() , idx.end() ) +1;
  F.allocate( nb_dof );
}

template<typename type>
template<typename T>
void
Builder<type>::transfer( const Field<type,T>& fx , Field<type,T>& fy ) const
{
  luma_implement;
}

} // luma
