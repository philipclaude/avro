#ifndef AVRO_MESH_DELAUNAY_H_
#define AVRO_MESH_DELAUNAY_H_

#include "mesh/points.h"
#include "common/types.h"

namespace avro
{

class Delaunay : public Points
{

public:
	Delaunay( const coord_t _dim ) : Points(_dim) {}
	Delaunay( const Points& points ) : Points(points) {}

	int bisector( const index_t i , const index_t j ) const;

	void seeds( const int b , index_t& i , index_t& j ) const;

	index_t closest( const real_t* x ) const;
};

} // avro

#endif
