// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef AVRO_MESH_DELAUNAY_H_
#define AVRO_MESH_DELAUNAY_H_

#include "mesh/vertices.h"
#include "common/types.h"

namespace avro
{

class Delaunay : public Vertices
{

public:
	Delaunay( const coord_t _dim ) : Vertices(_dim) {}
	Delaunay( const Vertices& vertices ) : Vertices(vertices) {}

	int bisector( const index_t i , const index_t j ) const;

	void seeds( const int b , index_t& i , index_t& j ) const;

	index_t closest( const real* x ) const;

};

}

#endif
