//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
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
  interpolate( const index_t* idx , index_t nv , const std::vector<real_t>& phi , type* u ) const;

};

} // avro

#endif
