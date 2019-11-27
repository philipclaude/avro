#ifndef LUNA_LIB_MESH_DOF_H_
#define LUNA_LIB_MESH_DOF_H_

#include "common/array.h"

namespace luna
{

template<typename type>
class DOF : public Array<type>
{
public:
  DOF( coord_t rank ) :
    Array<type>(ArrayLayout_Rectangular,rank)
  {}

  index_t rank() const { return Array<type>::rank_; }

  template<typename Shape>
  void
  interpolate( Shape& master , const index_t* idx , index_t nv , type* u )
  {
    u = type(0);
    for (index_t j=0;j<nv;j++)
    for (index_t d=0;d<rank();d++)
      u[d] += this->operator()(idx[j],d);
  }
};

} // luna

#endif
