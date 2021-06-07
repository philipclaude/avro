//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//

// no include guards because it should only be included once

#define USE_AVRO_PREDICATES 1
#define ALWAYS_EXACT 0
#define WITH_SIDE5 1

#if USE_AVRO_PREDICATES
	#include "numerics/predicates/side_exact.h"
	#include "numerics/predicates/side_filters.h"
#else
	#include "numerics/predicates/geo_side_exact.h"
	#include "numerics/predicates/geo_side_filters.h"
	using namespace GEO;
	using namespace PCK;

	namespace GEO
	{
		Sign ZERO = (Sign)avro::ZERO;
	}
#endif


#if USE_AVRO_PREDICATES
// exact definitions
#define side1_nd_exact_pck avro_side1_nd_exact_pck
#define side2_nd_exact_pck avro_side2_nd_exact_pck
#define side3_nd_exact_pck avro_side3_nd_exact_pck
#define side4_nd_exact_pck avro_side4_nd_exact_pck
#define side5_nd_exact_pck avro_side5_nd_exact_pck

// when ambient dimension = topological dimension
#define side4_3d_exact_pck avro_side4_3d_exact_pck

// filter definitions
#define side1_3d_filter avro_side1_3d_filter
#define side2_3d_filter avro_side2_3d_filter
#define side3_3d_filter avro_side3_3d_filter
#define side4_3d_filter avro_side4_3d_filter
#define side5_3d_filter avro_side5_3d_filter

#define side1_4d_filter avro_side1_4d_filter
#define side2_4d_filter avro_side2_4d_filter
#define side3_4d_filter avro_side3_4d_filter
#define side4_4d_filter avro_side4_4d_filter
#define side5_4d_filter avro_side5_4d_filter

#define side1_6d_filter avro_side1_6d_filter
#define side2_6d_filter avro_side2_6d_filter
#define side3_6d_filter avro_side3_6d_filter
#define side4_6d_filter avro_side4_6d_filter
#define side5_6d_filter avro_side5_6d_filter

#define side1_7d_filter avro_side1_7d_filter
#define side2_7d_filter avro_side2_7d_filter
#define side3_7d_filter avro_side3_7d_filter
#define side4_7d_filter avro_side4_7d_filter
#define side5_7d_filter avro_side5_7d_filter

#else
// exact definitions
#define side1_nd_exact_pck geo_side1_exact_SOS
#define side2_nd_exact_pck geo_side2_exact_SOS
#define side3_nd_exact_pck geo_side3_exact_SOS
#define side4_nd_exact_pck geo_side4_exact_SOS
#define side5_nd_exact_pck avro_side5_nd_exact_pck

// filter definitions
#define side1_3d_filter geo_side1_3d_filter
#define side2_3d_filter geo_side2_3d_filter
#define side3_3d_filter geo_side3_3d_filter
#define side4_3d_filter avro_side4_3d_filter
#define side5_3d_filter avro_side5_3d_filter

#define side1_4d_filter geo_side1_4d_filter
#define side2_4d_filter geo_side2_4d_filter
#define side3_4d_filter geo_side3_4d_filter
#define side4_4d_filter geo_side4_4d_filter
#define side5_4d_filter avro_side5_4d_filter

#define side1_6d_filter geo_side1_6d_filter
#define side2_6d_filter geo_side2_6d_filter
#define side3_6d_filter geo_side3_6d_filter
#define side4_6d_filter geo_side4_6d_filter
#define side5_6d_filter avro_side5_6d_filter

#define side1_7d_filter geo_side1_7d_filter
#define side2_7d_filter geo_side2_7d_filter
#define side3_7d_filter geo_side3_7d_filter
#define side4_7d_filter geo_side4_7d_filter
#define side5_7d_filter avro_side5_7d_filter
#endif

// these are not implemented in geogram
#define side1_5d_filter avro_side1_5d_filter
#define side2_5d_filter avro_side2_5d_filter
#define side3_5d_filter avro_side3_5d_filter
#define side4_5d_filter avro_side4_5d_filter
#define side5_5d_filter avro_side5_5d_filter

#define side1_8d_filter avro_side1_8d_filter
#define side2_8d_filter avro_side2_8d_filter
#define side3_8d_filter avro_side3_8d_filter
#define side4_8d_filter avro_side4_8d_filter
#define side5_8d_filter avro_side5_8d_filter

#define side1_9d_filter avro_side1_9d_filter
#define side2_9d_filter avro_side2_9d_filter
#define side3_9d_filter avro_side3_9d_filter
#define side4_9d_filter avro_side4_9d_filter
#define side5_9d_filter avro_side5_9d_filter

#define side1_10d_filter avro_side1_10d_filter
#define side2_10d_filter avro_side2_10d_filter
#define side3_10d_filter avro_side3_10d_filter
#define side4_10d_filter avro_side4_10d_filter
#define side5_10d_filter avro_side5_10d_filter

Sign side1_3d_SOS(const double* p0,const double* p1,const double* q0)
{
	Sign result = Sign(side1_3d_filter(p0,p1,q0));
	if(result==ZERO) {
		result = side1_nd_exact_pck(p0,p1,q0,3);
	}
	return result;
}

Sign side1_4d_SOS(const double* p0,const double* p1,const double* q0)
{
	Sign result = Sign(side1_4d_filter(p0,p1,q0));
	if(result==ZERO) {
		result = side1_nd_exact_pck(p0,p1,q0,4);
	}
	return result;
}

Sign side1_5d_SOS(const double* p0,const double* p1,const double* q0)
{
	Sign result = Sign(side1_5d_filter(p0,p1,q0));
	if(result==ZERO) {
		result = side1_nd_exact_pck(p0,p1,q0,5);
	}
	return result;
}

Sign side1_6d_SOS(const double* p0,const double* p1,const double* q0)
{
	Sign result = Sign(side1_6d_filter(p0,p1,q0));
	if(result==ZERO) {
		result = side1_nd_exact_pck(p0,p1,q0,6);
	}
	return result;
}

Sign side1_7d_SOS(const double* p0,const double* p1,const double* q0)
{
	Sign result = Sign(side1_7d_filter(p0,p1,q0));
	if(result==ZERO) {
		result = side1_nd_exact_pck(p0,p1,q0,7);
	}
	return result;
}

Sign side1_8d_SOS(const double* p0,const double* p1,const double* q0)
{
	Sign result = Sign(side1_8d_filter(p0,p1,q0));
	if(result==ZERO) {
		result = side1_nd_exact_pck(p0,p1,q0,8);
	}
	return result;
}

Sign side1_9d_SOS(const double* p0,const double* p1,const double* q0)
{
	Sign result = Sign(side1_9d_filter(p0,p1,q0));
	if(result==ZERO) {
		result = side1_nd_exact_pck(p0,p1,q0,9);
	}
	return result;
}

Sign side1_10d_SOS(const double* p0,const double* p1,const double* q0)
{
	Sign result = Sign(side1_10d_filter(p0,p1,q0));
	if(result==ZERO) {
		result = side1_nd_exact_pck(p0,p1,q0,10);
	}
	return result;
}

Sign side2_3d_SOS(const double* p0,const double* p1,const double* p2,const double* q0,const double* q1)
{
	Sign result = Sign(side2_3d_filter(p0,p1,p2,q0,q1));
	if(result==ZERO) {
		result = side2_nd_exact_pck(p0,p1,p2,q0,q1,3);
	}
	return result;
}

Sign side2_4d_SOS(const double* p0,const double* p1,const double* p2,const double* q0,const double* q1)
{
	Sign result = Sign(side2_4d_filter(p0,p1,p2,q0,q1));
	if(result==ZERO) {
		result = side2_nd_exact_pck(p0,p1,p2,q0,q1,4);
	}
	return result;
}

Sign side2_5d_SOS(const double* p0,const double* p1,const double* p2,const double* q0,const double* q1)
{
	Sign result = Sign(side2_5d_filter(p0,p1,p2,q0,q1));
	if(result==ZERO) {
		result = side2_nd_exact_pck(p0,p1,p2,q0,q1,5);
	}
	return result;
}

Sign side2_6d_SOS(const double* p0,const double* p1,const double* p2,const double* q0,const double* q1)
{
	Sign result = Sign(side2_6d_filter(p0,p1,p2,q0,q1));
	if(result==ZERO) {
		result = side2_nd_exact_pck(p0,p1,p2,q0,q1,6);
	}
	return result;
}

Sign side2_7d_SOS(const double* p0,const double* p1,const double* p2,const double* q0,const double* q1)
{
	Sign result = Sign(side2_7d_filter(p0,p1,p2,q0,q1));
	if(result==ZERO) {
		result = side2_nd_exact_pck(p0,p1,p2,q0,q1,7);
	}
	return result;
}

Sign side2_8d_SOS(const double* p0,const double* p1,const double* p2,const double* q0,const double* q1)
{
	Sign result = Sign(side2_8d_filter(p0,p1,p2,q0,q1));
	if(result==ZERO) {
		result = side2_nd_exact_pck(p0,p1,p2,q0,q1,8);
	}
	return result;
}

Sign side2_9d_SOS(const double* p0,const double* p1,const double* p2,const double* q0,const double* q1)
{
	Sign result = Sign(side2_9d_filter(p0,p1,p2,q0,q1));
	if(result==ZERO) {
		result = side2_nd_exact_pck(p0,p1,p2,q0,q1,9);
	}
	return result;
}

Sign side2_10d_SOS(const double* p0,const double* p1,const double* p2,const double* q0,const double* q1)
{
	Sign result = Sign(side2_10d_filter(p0,p1,p2,q0,q1));
	if(result==ZERO) {
		result = side2_nd_exact_pck(p0,p1,p2,q0,q1,10);
	}
	return result;
}

Sign side3_3d_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* q0,const double* q1,const double* q2)
{
	Sign result = Sign(side3_3d_filter(p0,p1,p2,p3,q0,q1,q2));
	if(result==ZERO) {
		result = side3_nd_exact_pck(p0,p1,p2,p3,q0,q1,q2,3);
	}
	return result;
}

Sign side3_4d_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* q0,const double* q1,const double* q2)
{
	Sign result = Sign(side3_4d_filter(p0,p1,p2,p3,q0,q1,q2));
	if(result==ZERO) {
		result = side3_nd_exact_pck(p0,p1,p2,p3,q0,q1,q2,4);
	}
	return result;
}

Sign side3_5d_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* q0,const double* q1,const double* q2)
{
	Sign result = Sign(side3_5d_filter(p0,p1,p2,p3,q0,q1,q2));
	if(result==ZERO) {
		result = side3_nd_exact_pck(p0,p1,p2,p3,q0,q1,q2,5);
	}
	return result;
}

Sign side3_6d_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* q0,const double* q1,const double* q2)
{
	Sign result = Sign(side3_6d_filter(p0,p1,p2,p3,q0,q1,q2));
	if(result==ZERO) {
		result = side3_nd_exact_pck(p0,p1,p2,p3,q0,q1,q2,6);
	}
	return result;
}

Sign side3_7d_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* q0,const double* q1,const double* q2)
{
	Sign result = Sign(side3_7d_filter(p0,p1,p2,p3,q0,q1,q2));
	if(result==ZERO) {
		result = side3_nd_exact_pck(p0,p1,p2,p3,q0,q1,q2,7);
	}
	return result;
}

Sign side3_8d_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* q0,const double* q1,const double* q2)
{
	Sign result = Sign(side3_8d_filter(p0,p1,p2,p3,q0,q1,q2));
	if(result==ZERO) {
		result = side3_nd_exact_pck(p0,p1,p2,p3,q0,q1,q2,8);
	}
	return result;
}

Sign side3_9d_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* q0,const double* q1,const double* q2)
{
	Sign result = Sign(side3_9d_filter(p0,p1,p2,p3,q0,q1,q2));
	if(result==ZERO) {
		result = side3_nd_exact_pck(p0,p1,p2,p3,q0,q1,q2,9);
	}
	return result;
}

Sign side3_10d_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* q0,const double* q1,const double* q2)
{
	Sign result = Sign(side3_10d_filter(p0,p1,p2,p3,q0,q1,q2));
	if(result==ZERO) {
		result = side3_nd_exact_pck(p0,p1,p2,p3,q0,q1,q2,10);
	}
	return result;
}

Sign side4_3d_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* p4,const double* q0,const double* q1,const double* q2,const double* q3)
{
	Sign result = Sign(side4_3d_filter(p0,p1,p2,p3,p4,q0,q1,q2,q3));
	if(result==ZERO) {
		//result = side4_3d_exact_pck(p0,p1,p2,p3,p4,true); // true for sos
		result = side4_nd_exact_pck(p0,p1,p2,p3,p4,q0,q1,q2,q3,3);
	}
	return result;
}

Sign side4_4d_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* p4,const double* q0,const double* q1,const double* q2,const double* q3)
{
	Sign result = Sign(side4_4d_filter(p0,p1,p2,p3,p4,q0,q1,q2,q3));
	if(result==ZERO) {
		result = side4_nd_exact_pck(p0,p1,p2,p3,p4,q0,q1,q2,q3,4);
	}
	return result;
}

Sign side4_5d_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* p4,const double* q0,const double* q1,const double* q2,const double* q3)
{
	Sign result = Sign(side4_5d_filter(p0,p1,p2,p3,p4,q0,q1,q2,q3));
	if(result==ZERO) {
		result = side4_nd_exact_pck(p0,p1,p2,p3,p4,q0,q1,q2,q3,5);
	}
	return result;
}

Sign side4_6d_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* p4,const double* q0,const double* q1,const double* q2,const double* q3)
{
	Sign result = Sign(side4_6d_filter(p0,p1,p2,p3,p4,q0,q1,q2,q3));
	if(result==ZERO) {
		result = side4_nd_exact_pck(p0,p1,p2,p3,p4,q0,q1,q2,q3,6);
	}
	return result;
}

Sign side4_7d_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* p4,const double* q0,const double* q1,const double* q2,const double* q3)
{
	Sign result = Sign(side4_7d_filter(p0,p1,p2,p3,p4,q0,q1,q2,q3));
	if(result==ZERO) {
		result = side4_nd_exact_pck(p0,p1,p2,p3,p4,q0,q1,q2,q3,7);
	}
	return result;
}

Sign side4_8d_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* p4,const double* q0,const double* q1,const double* q2,const double* q3)
{
	Sign result = Sign(side4_8d_filter(p0,p1,p2,p3,p4,q0,q1,q2,q3));
	if(result==ZERO) {
		result = side4_nd_exact_pck(p0,p1,p2,p3,p4,q0,q1,q2,q3,8);
	}
	return result;
}

Sign side4_9d_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* p4,const double* q0,const double* q1,const double* q2,const double* q3)
{
	Sign result = Sign(side4_9d_filter(p0,p1,p2,p3,p4,q0,q1,q2,q3));
	if(result==ZERO) {
		result = side4_nd_exact_pck(p0,p1,p2,p3,p4,q0,q1,q2,q3,9);
	}
	return result;
}

Sign side4_10d_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* p4,const double* q0,const double* q1,const double* q2,const double* q3)
{
	Sign result = Sign(side4_10d_filter(p0,p1,p2,p3,p4,q0,q1,q2,q3));
	if(result==ZERO) {
		result = side4_nd_exact_pck(p0,p1,p2,p3,p4,q0,q1,q2,q3,10);
	}
	return result;
}

#if WITH_SIDE5

Sign side5_4d_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* p4,const double* p5,const double* q0,const double* q1,const double* q2,const double* q3,const double* q4)
{
	Sign result = Sign(side5_4d_filter(p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4));
	if(result==ZERO) {
		result = side5_nd_exact_pck(p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4,4);
	}
	return result;
}

Sign side5_5d_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* p4,const double* p5,const double* q0,const double* q1,const double* q2,const double* q3,const double* q4)
{
	Sign result = Sign(side5_5d_filter(p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4));
	if(result==ZERO) {
		result = side5_nd_exact_pck(p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4,5);
	}
	return result;
}


Sign side5_6d_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* p4,const double* p5,const double* q0,const double* q1,const double* q2,const double* q3,const double* q4)
{
	Sign result = Sign(side5_6d_filter(p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4));
	if(result==ZERO) {
		result = side5_nd_exact_pck(p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4,6);
	}
	return result;
}

Sign side5_7d_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* p4,const double* p5,const double* q0,const double* q1,const double* q2,const double* q3,const double* q4)
{
	Sign result = Sign(side5_7d_filter(p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4));
	if(result==ZERO) {
		result = side5_nd_exact_pck(p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4,7);
	}
	return result;
}

Sign side5_8d_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* p4,const double* p5,const double* q0,const double* q1,const double* q2,const double* q3,const double* q4)
{
	Sign result = Sign(side5_8d_filter(p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4));
	if(result==ZERO) {
		result = side5_nd_exact_pck(p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4,8);
	}
	return result;
}

Sign side5_9d_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* p4,const double* p5,const double* q0,const double* q1,const double* q2,const double* q3,const double* q4)
{
	Sign result = Sign(side5_9d_filter(p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4));
	if(result==ZERO) {
		result = side5_nd_exact_pck(p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4,9);
	}
	return result;
}

Sign side5_10d_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* p4,const double* p5,const double* q0,const double* q1,const double* q2,const double* q3,const double* q4)
{
	Sign result = Sign(side5_10d_filter(p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4));
	if(result==ZERO) {
		result = side5_nd_exact_pck(p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4,10);
	}
	return result;
}

#endif

Sign side1_SOS(const double* p0,const double* p1,const double* q0, coord_index_t DIM)
{
	#if ALWAYS_EXACT
	return side1_nd_exact_pck(p0,p1,q0,DIM);
	#endif

	switch(DIM) {
	case 2:
		return side1_nd_exact_pck(p0,p1,q0,2);
	case 3:
		return side1_3d_SOS(p0,p1,q0);
	case 4:
		return side1_4d_SOS(p0,p1,q0);
	case 5:
		return side1_5d_SOS(p0,p1,q0);
	case 6:
		return side1_6d_SOS(p0,p1,q0);
	case 7:
		return side1_7d_SOS(p0,p1,q0);
	case 8:
		return side1_8d_SOS(p0,p1,q0);
	case 9:
		return side1_9d_SOS(p0,p1,q0);
	case 10:
		return side1_10d_SOS(p0,p1,q0);
	}
	geo_assert_not_reached;
	return GEO::ZERO; // suppress compiler warning
}

Sign side2_SOS(const double* p0,const double* p1,const double* p2,const double* q0,const double* q1, coord_index_t DIM)
{
	#if ALWAYS_EXACT
	return side2_nd_exact_pck(p0,p1,p2,q0,q1,DIM);
	#endif

	switch(DIM) {
	case 2:
		return side2_nd_exact_pck(p0,p1,p2,q0,q1,2);
	case 3:
		return side2_3d_SOS(p0,p1,p2,q0,q1);
	case 4:
		return side2_4d_SOS(p0,p1,p2,q0,q1);
	case 5:
		return side2_5d_SOS(p0,p1,p2,q0,q1);
	case 6:
		return side2_6d_SOS(p0,p1,p2,q0,q1);
	case 7:
		return side2_7d_SOS(p0,p1,p2,q0,q1);
	case 8:
		return side2_8d_SOS(p0,p1,p2,q0,q1);
	case 9:
		return side2_9d_SOS(p0,p1,p2,q0,q1);
	case 10:
		return side2_10d_SOS(p0,p1,p2,q0,q1);
	}
	geo_assert_not_reached;
	return GEO::ZERO; // suppress compiler warning
}

Sign side3_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* q0,const double* q1,const double* q2, coord_index_t DIM)
{
	#if ALWAYS_EXACT
	return side3_nd_exact_pck(p0,p1,p2,p3,q0,q1,q2,DIM);
	#endif

	switch(DIM) {
	case 2:
		return side3_nd_exact_pck(p0,p1,p2,p3,q0,q1,q2,2);
	case 3:
		return side3_3d_SOS(p0,p1,p2,p3,q0,q1,q2);
	case 4:
		return side3_4d_SOS(p0,p1,p2,p3,q0,q1,q2);
	case 5:
		return side3_5d_SOS(p0,p1,p2,p3,q0,q1,q2);
	case 6:
		return side3_6d_SOS(p0,p1,p2,p3,q0,q1,q2);
	case 7:
		return side3_7d_SOS(p0,p1,p2,p3,q0,q1,q2);
	case 8:
		return side3_8d_SOS(p0,p1,p2,p3,q0,q1,q2);
	case 9:
		return side3_9d_SOS(p0,p1,p2,p3,q0,q1,q2);
	case 10:
		return side3_10d_SOS(p0,p1,p2,p3,q0,q1,q2);
	}
	geo_assert_not_reached;
	return GEO::ZERO; // suppress compiler warning
}

Sign side4_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* p4,const double* q0,const double* q1,const double* q2,const double* q3, coord_index_t DIM)
{
	#if ALWAYS_EXACT
	return side4_nd_exact_pck(p0,p1,p2,p3,p4,q0,q1,q2,q3,DIM);
	#endif

	switch(DIM) {
	case 3:
		return side4_3d_SOS(p0,p1,p2,p3,p4,q0,q1,q2,q3);
	case 4:
		return side4_4d_SOS(p0,p1,p2,p3,p4,q0,q1,q2,q3);
	case 5:
		return side4_5d_SOS(p0,p1,p2,p3,p4,q0,q1,q2,q3);
	case 6:
		return side4_6d_SOS(p0,p1,p2,p3,p4,q0,q1,q2,q3);
	case 7:
		return side4_7d_SOS(p0,p1,p2,p3,p4,q0,q1,q2,q3);
	case 8:
		return side4_8d_SOS(p0,p1,p2,p3,p4,q0,q1,q2,q3);
	case 9:
		return side4_9d_SOS(p0,p1,p2,p3,p4,q0,q1,q2,q3);
	case 10:
		return side4_10d_SOS(p0,p1,p2,p3,p4,q0,q1,q2,q3);
	}
	geo_assert_not_reached;
	return GEO::ZERO; // suppress compiler warning
}

Sign side5_SOS(const double* p0,const double* p1,const double* p2,const double* p3,const double* p4,const double* p5,const double* q0,const double* q1,const double* q2,const double* q3,const double* q4, coord_index_t DIM)
{
	#if ALWAYS_EXACT
	return side5_nd_exact_pck(p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4,DIM);
	#endif

	switch(DIM) {
	case 4:
		return side5_4d_SOS(p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4);
	case 5:
		return side5_5d_SOS(p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4);
	case 6:
		return side5_6d_SOS(p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4);
	case 7:
		return side5_7d_SOS(p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4);
	case 8:
		return side5_8d_SOS(p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4);
	case 9:
		return side5_9d_SOS(p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4);
	case 10:
		return side5_10d_SOS(p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4);
	}
	geo_assert_not_reached;
	return GEO::ZERO; // suppress compiler warning
}
