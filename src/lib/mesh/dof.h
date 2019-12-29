#ifndef luma_LIB_MESH_DOF_H_
#define luma_LIB_MESH_DOF_H_

#include "common/table.h"

namespace luma
{

template<typename type>
class DOF : public Table<type>
{
public:
  DOF( coord_t rank ) :
    Table<type>(TableLayout_Rectangular,rank)
  {}

  index_t rank() const { return Table<type>::rank_; }

  void
  interpolate( const std::vector<const type*>& uk , const std::vector<real_t>& phi , type* u ) const
  {
    luma_assert( uk.size() == phi.size() );
    for (index_t j=0;j<rank();j++)
      u[j] = type(0);
    for (index_t j=0;j<rank();j++)
    {
      for (index_t i=0;i<uk.size();i++)
        u[j] += uk[i][j]*phi[i];
    }
  }

  void
  interpolate( const index_t* idx , index_t nv , const std::vector<real_t>& phi , type* u ) const
  {
    luma_assert( nv == phi.size() );
    for (index_t j=0;j<rank();j++)
      u[j] = type(0);
    for (index_t j=0;j<rank();j++)
    {
      for (index_t i=0;i<phi.size();i++)
        u[j] += (*this)( idx[i] , 0 )*phi[i];
    }
  }

};

} // luma

#endif
