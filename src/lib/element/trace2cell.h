//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef AVRO_LIB_ELEMENT_TRACE2CELL_H_
#define AVRO_LIB_ELEMENT_TRACE2CELL_H_

#include "avro_types.h"
#include "common/error.h"

#include "element/simplex.h"

namespace avro
{

class CanonicalTraceToCell
{
public:
  CanonicalTraceToCell( index_t _elem , index_t _trace , int _orient ) :
        elem(_elem),
        trace(_trace),
        orient(_orient)
  {}

  index_t elem;
  index_t trace;
  int orient;
};

class TraceToCellRefCoord
{
public:
  TraceToCellRefCoord( const Simplex& trace ) :
    trace_(trace)
  {}

  void eval( const CanonicalTraceToCell& canonical , const real_t* ut , real_t* uc ) const;

private:
  const Simplex& trace_;

};

} // avro

#endif
