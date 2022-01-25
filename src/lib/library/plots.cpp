//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "library/plots.h"

#include "element/simplex.h"

#include "mesh/points.h"

namespace avro
{

namespace library
{

template<typename type>
Plot<type>::Plot( Points& points ) :
  Topology<type>(points,0)
{
  for (index_t k=0;k<points.nb();k++)
    Table<index_t>::add( &k , 1 );
}

template class Plot<Simplex>;

} // library

} // avro
