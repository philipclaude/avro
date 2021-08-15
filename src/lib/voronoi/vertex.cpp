//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "common/set.h"

#include "element/simplex.h"
#include "element/polytope.h"

#include "mesh/points.h"

#include "numerics/predicates.h"
#include "voronoi/vertex.h"

#include <triangle/predicates.h>

#include <cmath>

namespace avro
{

namespace voronoi
{

Vertex::Vertex( const coord_t _dim ) :
	dim_(_dim),
	number_(0),
	x_(dim_,0.0)
{}

Vertex::Vertex( const coord_t _dim , const coord_t _number ) :
	dim_(_dim),
	number_(_number),
	x_(dim_,0.0)
{}

Vertex::Vertex( const std::vector<real_t>& _x , const coord_t _number ) :
	dim_( _x.size() ) ,
	number_(_number) ,
	x_(_x)
{}

Vertex::Vertex( const real_t* _x , const coord_t _dim ) :
	dim_(_dim),
	number_(0),
	x_(_x,_x+dim_)
{}

Vertex::Vertex( const Vertex& v0 ) :
	dim_( v0.dim() ),
	number_(v0.number()),
	x_( v0.x() ),
	bisector_( v0.bisectors() ) ,
	simplex_( v0.simplices() ) ,
	site_( v0.sites() )
{}

void
Vertex::add_site( const real_t* zj ) {
	site_.push_back(zj);
}

void
Vertex::intersect_symbolic( const Vertex* v , const Vertex* v0 ) {
	intersect_bisectors(v,v0);
	intersect_simplices(v,v0);
}

void
Vertex::intersect_bisectors( const Vertex* v0 , const Vertex* v1 ) {
	Set::intersection( v0->bisectors() , v1->bisectors() , bisector_ );
	avro_assert( bisector_.size()==index_t(number_-1) );
}

void
Vertex::intersect_simplices( const Vertex* v0 , const Vertex* v1 ) {
	#if 0
	Set::intersection( v0->simplices() , v1->simplices() , simplex_ );
	uniquify(simplex_);
	#else
	simplex_.resize( v0->nb_simplices() +v1->nb_simplices() );
	index_t j = 0;
	for (index_t k = 0; k < v0->nb_simplices(); k++, j++)
		simplex_[k] = v0->simplex(k);
	for (index_t k = 0; k < v1->nb_simplices(); k++, j++)
		simplex_[j] = v1->simplex(k);
	uniquify(simplex_);
	#endif
}

void
Vertex::set_sites( const Points& delaunay , const std::map<int,Bisector>& B ) {

	index_t Pi_i, Pi_j;
	site_.clear();
	for (index_t k = 0; k < bisector_.size(); k++) {
		if (bisector_[k] < 0) continue; // skip mesh facets

		// get the actual seed indices in the delaunay container
		std::map<int,Bisector>::const_iterator it = B.find(bisector_[k]);
		avro_assert( it != B.end() );

		Pi_i = it->second.p0;
		Pi_j = it->second.p1;
		if ( delaunay[Pi_j] == z0_ )
			site_.push_back( delaunay[Pi_i] );
		else
			site_.push_back( delaunay[Pi_j] );
	}
	uniquify(site_);
}

void
Vertex::intersect_geometric( const real_t* q1 , const real_t* q2 , const real_t* p1 , const real_t *p2 ) {

	double d = 0., l1 = 0., l2 = 0.;
  for (coord_t c = 0; c < dim_; c++) {
    double n = p1[c] -p2[c];
    d -= n*(p2[c] +p1[c]);
    l1 += q2[c]*n;
    l2 += q1[c]*n;
  }
  d = .5*d;
  l1 = ::fabs(l1 +d);
  l2 = ::fabs(l2 +d);
  double l12 = l1 +l2;
  if (l12 > 1e-30) {
    l1 /= l12;
    l2 /= l12;
  }
  else {
    l1 = .5;
    l2 = .5;
  }

	// set the coordinates
  for (coord_t c = 0; c < dim_; c++)
    x_[c] = l1*q1[c] +l2*q2[c];
}

GEO::Sign
Vertex::side_inexact( const real_t *zi , const real_t *zj ) {

  double di = 0.,dj = 0.;
  for (coord_t d = 0; d < dim_; d++) {
    di += ( x_[d] -zi[d] )*( x_[d] -zi[d] );
    dj += ( x_[d] -zj[d] )*( x_[d] -zj[d] );
  }
  if (di > dj) return GEO::NEGATIVE;
  return GEO::POSITIVE;
}

GEO::Sign
Vertex::side(const real_t *zi , const real_t *zj , const bool exact )
{
  // fast version
	GEO::Sign result;
  if (!exact) {
    return side_inexact( zi , zj );
  }

	if (simplex_.size() != site_.size()+1)
		print("v",true);

  avro_assert_msg( simplex_.size() == site_.size()+1 ,
	 								 "|s| = %lu, nb_sites = %lu",simplex_.size(),site_.size() );

  switch( simplex_.size() ) {
    case 1:
      // this vertex did not result from an intersection, just determine
      // which side it lies on normally
      result = GEO::PCK::side1_SOS(zi, zj, simplex_[0], dim());
      avro_assert( result != GEO::ZERO );
      break;
    case 2:
      // this vertex is the intersection of a bisector with an edge (s0-s1) of the topology
      result = GEO::PCK::side2_SOS(zi, site_[0], zj, simplex_[0], simplex_[1], dim_);
      avro_assert( result != GEO::ZERO );
      break;
    case 3:
      // this vertex is the intersection of two bisectors with a triangle (p0-p1-p2) of the topology
			try {
      	result = GEO::PCK::side3_SOS(zi, site_[0], site_[1], zj, simplex_[0], simplex_[1], simplex_[2], dim_);
			}
			catch(...) {
				print("v",true);
				avro_implement;
			}
			if (result == GEO::ZERO) {
				for (index_t k = 0; k < site_.size(); k++) {
					std::vector<real_t> simplex_coords(simplex_[k],simplex_[k]+dim_);
					print_inline( simplex_coords , "simplex " + stringify(k ) );
				}

				std::vector<real_t> zi_coord( zi , zi+dim_ );
				std::vector<real_t> zj_coord( zj , zj+dim_ );
				print_inline( zi_coord , "zi" );
				print_inline( zj_coord , "zj" );

			}
      avro_assert( result != GEO::ZERO );
      break;
    case 4:
      // this vertex is the intersection of three bisectors with a tetrahedron (p0-p1-p2-p3) of the topology
      result = GEO::PCK::side4_SOS(zi, site_[0], site_[1], site_[2], zj, simplex_[0], simplex_[1], simplex_[2], simplex_[3], dim_);
			if (result == GEO::ZERO) {
				for (index_t k = 0; k < simplex_.size(); k++) {
					std::vector<real_t> simplex_coords(simplex_[k],simplex_[k]+dim_);
					print_inline( simplex_coords , "simplex " + stringify(k ) );
				}
				for (index_t k = 0; k < site_.size(); k++) {
					std::vector<real_t> site_coords(simplex_[k],simplex_[k]+dim_);
					print_inline( site_coords , "site " + stringify(k ) );
				}
				std::vector<real_t> zi_coord( zi , zi+dim_ );
				std::vector<real_t> zj_coord( zj , zj+dim_ );
				print_inline( zi_coord , "zi" );
				print_inline( zj_coord , "zj" );
				print("v",true);
				//printf("zero = %d, geo::zero = %d, equal? %s\n",ZERO,GEO::ZERO,ZERO==GEO::ZERO?"true":"false");
			}
      avro_assert( result != GEO::ZERO );
      break;
    case 5:
      // this vertex is the intersection of four bisectors with a pentatope (p0-p1-p2-p3-p4) of the topology
      result = GEO::PCK::side5_SOS(zi, site_[0], site_[1], site_[2], site_[3], zj, simplex_[0], simplex_[1], simplex_[2], simplex_[3], simplex_[4], dim_);
			if (result == GEO::ZERO) return result;
      //avro_assert( result!=GEO::ZERO );
      break;
    default:
      printf("more predicates needed for %d-simplices!",(int)simplex_.size());
      avro_implement;
  }

  avro_assert( result != GEO::ZERO );

  return result;
}

void
Vertex::print( const std::string& pre , const bool symbolic ) const {

	// print the coordinates
	printf("%s: (",pre.c_str());
	for (coord_t d = 0; d < dim_; d++)
		printf(" %3.2e",x_[d]);
	printf(" )");

	// option to print symbolic information
	if (symbolic) {

		// print the bisectors
		printf(" - [");
		for (index_t k = 0; k < bisector_.size(); k++)
			printf(" b%d ",bisector_[k]);
		printf("]");

		// print the simplex vertex information
		printf(" - [");
		for (index_t k = 0; k < simplex_.size(); k++)
			printf(" s%p ",(void*)simplex_[k]);
		printf("]");

		// print the sites
		printf(" - [");
		for (index_t k = 0; k < site_.size(); k++)
			printf(" z%p ",(void*)site_[k]);
		printf("]");
	}
	printf("\n");
}

} // voronoi

} // avro
