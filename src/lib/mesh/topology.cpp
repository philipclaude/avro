// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "mesh/topology.h"
#include "mesh/vertices.h"

namespace ursa
{

void
TopologyBase::copy( TopologyBase& topology1 )
{
  dummy_    = topology1.dummy();
  number_   = topology1.number();
  name_     = topology1.name();
  elements_ = topology1.elements();
  first_    = topology1.first();
  last_     = topology1.last();
  vertices_ = topology1.vertices();
}

} // avro
