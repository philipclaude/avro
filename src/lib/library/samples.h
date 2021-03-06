//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_LIBRARY_SAMPLES_H_
#define avro_LIB_LIBRARY_SAMPLES_H_

#include "mesh/topology.h"
#include "mesh/points.h"

namespace avro
{

namespace library
{

class TwoTriangles : public Topology<Simplex>
{
public:
  TwoTriangles();

  Topology<Simplex>& edges() { return edges_; }

private:
  Points points_;
  Topology<Simplex> edges_;
};

class RegularPolygon : public Topology<Polytope>
{
public:
  RegularPolygon( index_t nb_side );
private:
  Points points_;  
};

} // library

} // avro

#endif
