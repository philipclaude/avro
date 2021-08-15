//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef AVRO_MESH_DELAUNAY_H_
#define AVRO_MESH_DELAUNAY_H_

#include "mesh/points.h"
#include "avro_types.h"

namespace avro
{

#if 0
class Delaunay : public Points
{

public:
	Delaunay( const coord_t _dim ) : Points(_dim) {}
	Delaunay( const Points& points ) : Points(points) {}

	int bisector( const index_t i , const index_t j ) const;

	void seeds( const int b , index_t& i , index_t& j ) const;

	index_t closest( const real_t* x ) const;
};

#endif

} // avro

#endif
