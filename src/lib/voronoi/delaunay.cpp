#include "common/error.h"

#include "mesh/points.h"

#include "voronoi/delaunay.h"

#include "numerics/geometry.h"

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
Delaunay::closest( const real_t* x ) const
{
	avro_assert( nb()>0 );
	int k0 = 0;
	real_t d,dmin = numerics::distance2( x , operator[](0) , dim_ );
	for (index_t k=1;k<nb();k++)
	{
		d = numerics::distance2( x , operator[](k) , dim_ );
		if (d<dmin)
		{
			dmin = d;
			k0   = k;
		}
	}
	return index_t(k0);
}

}
