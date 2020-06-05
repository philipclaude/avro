//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_NUMERICS_COORDINATE_H_
#define avro_LIB_NUMERICS_COORDINATE_H_

#include "common/types.h"

#include <vector>

namespace avro
{

namespace numerics
{

class Coordinate : public std::vector<real_t>
{
public:
  using std::vector<real_t>::data;
  using std::vector<real_t>::operator[];

  Coordinate( const coord_t dim ) :
    std::vector<real_t>(dim)
  {}

  Coordinate(real_t* x , const coord_t dim ) :
    std::vector<real_t>(x,x+dim)
  {}

  coord_t dim() const { return size(); }

private:
  using std::vector<real_t>::size;

};

} // numerics

} // avro

#endif
