#ifndef URSA_LIB_NUMERICS_GEOMETRY_H_
#define URSA_LIB_NUMERICS_GEOMETRY_H_

#include "common/types.h"

#include <vector>
#include <math.h>

namespace ursa
{

class Points;

namespace numerics
{

real_t distanceLinf( const real_t* x , const real_t* y , const coord_t dim );
real_t distanceLinf2( const real_t* x , const real_t* y , const coord_t dim );
real_t distance2( const real_t* x , const real_t* y , const coord_t dim );
real_t distance( const real_t* x , const real_t* y , const coord_t dim );
real_t volume( const std::vector<real_t*>& x , const coord_t dim );
real_t signedVolume( const std::vector<real_t*>& x , const coord_t dim );
real_t volume( const std::vector<const real_t*>& x , const coord_t dim );
real_t volume2( const std::vector<real_t*>& x );
real_t volume3( const std::vector<real_t*>& x );
real_t volume4( const std::vector<real_t*>& x );
real_t orient4d( const real_t* , const real_t* , const real_t* , const real_t* , const real_t* );

void centroid( const index_t* v0 , const index_t nv , const Points& v , std::vector<real_t>& xc );
void vector( const real_t* v0 , const real_t* v1 , const coord_t dim , real_t* v );

void normal( const std::vector<real_t*>& x , real_t* n , const coord_t dim );
void orthogonal( const std::vector<real_t*>& x , real_t* n , const coord_t dim );
void orientNormal( const std::vector<real_t*>& x , real_t* n , const coord_t dim );
void normalize( real_t* x , const coord_t dim );
void axpb( const real_t a , const real_t* x , const real_t* b , const coord_t dim , real_t* y );

void barycentric( const real_t* p , const std::vector<const real_t*>& x , const coord_t dim , std::vector<real_t>& alpha );
void barycentric_signed( real_t* p , const std::vector<real_t*>& x , const coord_t dim , std::vector<real_t>& alpha );

} // numerics

} // ursa

#endif
