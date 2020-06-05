//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_LIBRARY_EGADS_H_
#define avro_LIB_LIBRARY_EGADS_H_

#include "geometry/egads/body.h"

namespace avro
{

namespace EGADS
{

class Cube : public Body
{
public:
  Cube( const Context* context , coord_t dim );
  Cube( const Context* context , const std::vector<real_t>& lens , const real_t* x0=nullptr );

private:
  ego object_;
};

} // EGADS

} // avro

#endif
