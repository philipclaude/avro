#ifndef avro_LIB_MESH_DOF_H_
#define avro_LIB_MESH_DOF_H_

#include "common/table.h"

namespace avro
{

template<typename type>
class DOF : public Table<type>
{
public:
  DOF( coord_t rank ) :
    Table<type>(TableLayout_Rectangular,rank)
  {}

  index_t rank() const { return Table<type>::rank_; }

  bool
  interpolate( const std::vector<const type*>& uk , const std::vector<real_t>& phi , type* u ) const
  {
    avro_assert( uk.size() == phi.size() );
    for (index_t j=0;j<rank();j++)
      u[j] = type(0);
    for (index_t j=0;j<rank();j++)
    {
      for (index_t i=0;i<uk.size();i++)
        u[j] += uk[i][j]*phi[i];
    }
    return true;
  }

  bool
  interpolate( const index_t* idx , index_t nv , const std::vector<real_t>& phi , type* u ) const
  {
    avro_assert_msg( nv == phi.size() , "nv = %lu, |phi| = %lu" , nv , phi.size() );
    for (index_t j=0;j<rank();j++)
      u[j] = type(0);
    for (index_t j=0;j<rank();j++)
    {
      for (index_t i=0;i<phi.size();i++)
        u[j] += (*this)( idx[i] , 0 )*phi[i];
    }
    return true;
  }

};

} // avro

#endif
