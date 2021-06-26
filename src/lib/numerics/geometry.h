//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_NUMERICS_GEOMETRY_H_
#define avro_LIB_NUMERICS_GEOMETRY_H_

#include "types.h"

#include <vector>
#include <cmath>

namespace avro
{

class Points;

namespace numerics
{

real_t distanceLinf( const real_t* x , const real_t* y , const coord_t dim );
real_t distanceLinf2( const real_t* x , const real_t* y , const coord_t dim );
real_t distance2( const real_t* x , const real_t* y , const coord_t dim );
real_t distance( const real_t* x , const real_t* y , const coord_t dim );
real_t volume( const std::vector<const real_t*>& x , const coord_t dim );
real_t simplex_volume( const std::vector<const real_t*>& x , const coord_t dim );
real_t volume_nd( const std::vector<const real_t*>& x , const coord_t dim );
real_t volume2( const std::vector<const real_t*>& x );
real_t volume3( const std::vector<const real_t*>& x );
real_t volume4( const std::vector<const real_t*>& x );
real_t orient4d( const real_t* , const real_t* , const real_t* , const real_t* , const real_t* );

void centroid( const index_t* v0 , const index_t nv , const Points& v , std::vector<real_t>& xc );
void vector( const real_t* v0 , const real_t* v1 , const coord_t dim , real_t* v );

void normal( const std::vector<real_t*>& x , real_t* n , const coord_t dim );
void orthogonal( const std::vector<real_t*>& x , real_t* n , const coord_t dim );
void orientNormal( const std::vector<real_t*>& x , real_t* n , const coord_t dim );
void normalize( real_t* x , const coord_t dim );
void axpb( const real_t a , const real_t* x , const real_t* b , const coord_t dim , real_t* y );

void barycentric( const real_t* p , const std::vector<const real_t*>& x , const coord_t dim , std::vector<real_t>& alpha );
void barycentric_signed( real_t* p , const std::vector<const real_t*>& x , const coord_t dim , std::vector<real_t>& alpha );

} // numerics

} // avro

#endif
