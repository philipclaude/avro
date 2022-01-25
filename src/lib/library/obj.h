//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_LIBRARY_OBJ_H_
#define avro_LIB_LIBRARY_OBJ_H_

#include "mesh/topology.h"
#include "mesh/points.h"

#include <string>

namespace avro
{

namespace library
{

class objFile : public Topology<Simplex>
{
public:
  objFile( const std::string& filename );

  void read();

private:
  Points points_;
  std::string filename_;
};

} // library

} // avro

#endif
