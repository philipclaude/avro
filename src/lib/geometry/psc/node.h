//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_GEOMETRY_PSC_NODE_H_
#define avro_LIB_GEOMETRY_PSC_NODE_H_

#include "geometry/psc/object.h"

namespace avro
{

namespace PSC
{

class Node : public Object
{
public:
  Node( Body* body , real_t* data ) :
    Object(0,body->dim()),
    x_(data,data+dim_)
  {}

  real_t operator()(coord_t d) const
  {
    avro_assert( d<dim_ );
    return x_[d];
  }

  real_t& operator()(coord_t d)
  {
    avro_assert( d<dim_ );
    return x_[d];
  }

  const real_t* x() const { return x_.data(); }

  void inverse( std::vector<real_t>& x , std::vector<real_t>& u ) const;
  void inverse_guess( std::vector<real_t>& x , std::vector<real_t>& u ) const
  { inverse(x,u); }
  void evaluate( const std::vector<real_t>& u , std::vector<real_t>& x ) const;


private:
  std::vector<real_t> x_;
};

} // PSC

} // avro

#endif
