// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "common/error.h"

#include "mesh/geometrics.h"
#include "mesh/vertices.h"

#include "mesh/delaunay/delaunay.h"

namespace avro
{

int
Delaunay::bisector( const index_t i , const index_t j ) const
{
	// ensure the same bisector label is used for both sides
	if (i<j)
		return i +j*nb();
	return j +i*nb();
}

void
Delaunay::seeds( const int b , index_t& i , index_t& j ) const
{
	j = (index_t) b/nb();
	i = b -j*nb();
}

index_t
Delaunay::closest( const real* x ) const
{
	avro_assert( nb()>0 );
	int k0 = 0;
	real d,dmin = geometrics::distance2( x , operator[](0) , dim_ );
	for (index_t k=1;k<nb();k++)
	{
		d = geometrics::distance2( x , operator[](k) , dim_ );
		if (d<dmin)
		{
			dmin = d;
			k0   = k;
		}
	}
	return index_t(k0);
}

}
