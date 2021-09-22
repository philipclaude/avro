//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_LIBRARY_PLOTS_H_
#define avro_LIB_LIBRARY_PLOTS_H_

#include "mesh/topology.h"

namespace avro
{

class Points;

namespace library
{

template<typename type>
class Plot : public Topology<type>
{
public:
  Plot( Points& points );
  Plot( Topology<type>& topology );
};

} // library

} // avro

#endif
