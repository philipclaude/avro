//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_VORONOI_CDT_H_
#define avro_LIB_VORONOI_CDT_H_

#include "element/simplex.h"

#include "mesh/topology.h"

namespace avro
{

class PSC;

class ConstrainedDelaunayTriangulation : public Topology<Simplex>
{
public:
  ConstrainedDelaunayTriangulation( const Topology<Simplex>& psc );
  ConstrainedDelaunayTriangulation( const PSC& psc ); // will need to be broken into points, segments, triangles, etc.

  void compute();

  void recover_segments();
  void recover_facets();

private:
  const Topology<Simplex>& psc_;
};

} // avro

#endif
