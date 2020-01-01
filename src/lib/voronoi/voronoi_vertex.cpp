#include "common/set.h"

#include "master/simplex.h"
#include "master/polytope.h"

#include "mesh/points.h"

#include "numerics/predicates.h"

#include "voronoi/delaunay.h"
#include "voronoi/voronoi.h"

#include <cmath>

namespace avro
{

namespace delaunay
{

Vertex::Vertex( const coord_t _dim ) :
	dim_(_dim) , number_(0)
{
	init();
}

Vertex::Vertex( const coord_t _dim , const coord_t _number ) :
	dim_(_dim) , number_(_number)
{
	init();
}

Vertex::Vertex( const std::vector<real_t> _x , const coord_t _number ) :
	dim_( _x.size() ) , number_(_number) , x_(_x)
{
}

Vertex::Vertex( const real_t* _x , const coord_t _dim ) :
	dim_(_dim)
{
	init();
	for (coord_t d=0;d<dim_;d++)
		x_[d] = _x[d];
}

Vertex::Vertex( const Vertex& v0 ) :
	dim_( v0.dim() ) , x_( v0.x() ) , bisector_( v0.bisectors() ) ,
	topology_( v0.topologies() ) , simplex_( v0.simplices() ) , site_( v0.sites() )
{
}

void
Vertex::init()
{
	x_.resize(dim_);
	for (coord_t d=0;d<dim_;d++)
		x_[d] = 0.;
}

void
Vertex::intersectSymbolic( const Vertex* v , const Vertex* v0 ,
	const Delaunay& delaunay )
{
	intersectBisectors(v,v0);
	intersectMeshes(v,v0);
	intersectSimplices(v,v0);
	setSites(delaunay);
}

void
Vertex::intersectBisectors( const Vertex* v0 , const Vertex* v1 )
{
	Set::intersection( v0->bisectors() , v1->bisectors() , bisector_ );
	avro_assert( bisector_.size()==index_t(number_-1) );
}

void
Vertex::intersectMeshes( const Vertex* v0 , const Vertex* v1 )
{
	Set::intersection( v0->topologies() , v1->topologies() , topology_ );
}

void
Vertex::intersectSimplices( const Vertex* v0 , const Vertex* v1 )
{
	//Set::intersection( v0->simplices() , v1->simplices() , simplex_ );
	simplex_.resize( v0->nb_simplices() +v1->nb_simplices() );
	index_t j=0;
	for (index_t k=0;k<v0->nb_simplices();k++,j++)
		simplex_[k] = v0->simplex(k);
	for (index_t k=0;k<v1->nb_simplices();k++,j++)
		simplex_[j] = v1->simplex(k);
	uniquify(simplex_);
}

void
Vertex::setSites( const Delaunay& delaunay )
{
	index_t Pi_i, Pi_j;

	site_.clear();

	for (index_t k=0;k<bisector_.size();k++)
	{
		if (bisector_[k]<0) continue; // skip mesh facets

		// get the actual seed indices in the delaunay container
		delaunay.seeds( bisector_[k] , Pi_i , Pi_j );

		if ( delaunay[Pi_j]==z0_ )
			site_.push_back( delaunay[Pi_i] );
		else
			site_.push_back( delaunay[Pi_j] );
	}
	uniquify(site_);
}

void
Vertex::intersectGeometric( const real_t* q1 , const real_t* q2 , const real_t* p1 , const real_t *p2 )
{
	double d = 0., l1 = 0., l2 = 0.;
  for (coord_t c=0;c<dim_;c++)
  {
    double n = p1[c] -p2[c];
    d -= n*(p2[c] +p1[c]);
    l1 += q2[c]*n;
    l2 += q1[c]*n;
  }
  d = .5*d;
  l1 = ::fabs(l1 +d);
  l2 = ::fabs(l2 +d);
  double l12 = l1 +l2;
  if (l12>1e-30)
  {
    l1 /= l12;
    l2 /= l12;
  }
  else
  {
    l1 = .5;
    l2 = .5;
  }

	// set the coordinates
  for (coord_t c=0;c<dim_;c++)
    x_[c] = l1*q1[c] +l2*q2[c];
}

GEO::Sign
Vertex::sideFast( const real_t *zi , const real_t *zj )
{
  double di=0.,dj=0.;
  for (coord_t d=0;d<dim_;d++)
  {
    di += ( x_[d] -zi[d] )*( x_[d] -zi[d] );
    dj += ( x_[d] -zj[d] )*( x_[d] -zj[d] );
  }
  if (di>dj) return GEO::NEGATIVE;
  return GEO::POSITIVE;
}

GEO::Sign
Vertex::side(const real_t *zi , const real_t *zj , const bool exact )
{
  GEO::Sign result;

  // fast version
  if (!exact)
  {
    return sideFast( zi , zj );
  }

	if (simplex_.size()!=site_.size()+1)
	{
		print("verror",true);
	}

  avro_assert_msg( simplex_.size()==site_.size()+1 ,
	 								 "|s| = %lu, nb_sites = %lu",simplex_.size(),site_.size() );

  switch( simplex_.size() )
  {
    case 1:
      // this vertex did not result from an intersection, just determine
      // which side it lies on normally
      result = GEO::PCK::side1_SOS(zi,zj,simplex_[0],dim());
      avro_assert( result!=GEO::ZERO );
      break;
    case 2:
      // this vertex is the intersection of a bisector with an edge (s0-s1) of the topology
      result = GEO::PCK::side2_SOS(zi,site_[0],zj,simplex_[0],simplex_[1],dim_);
      avro_assert( result!=GEO::ZERO );
      break;
    case 3:
      // this vertex is the intersection of two bisectors with a triangle (p0-p1-p2) of the topology
      result = GEO::PCK::side3_SOS(zi,site_[0],site_[1],zj,simplex_[0],simplex_[1],simplex_[2],dim_);
      avro_assert( result!=GEO::ZERO );
      break;
    case 4:
      // this vertex is the intersection of three bisectors with a tetrahedron (p0-p1-p2-p3) of the topology
      result = GEO::PCK::side4_SOS(zi,site_[0],site_[1],site_[2],zj,simplex_[0],simplex_[1],simplex_[2],simplex_[3],dim_);
      avro_assert( result!=GEO::ZERO );
      break;
    case 5:
      // this vertex is the intersection of four bisectors with a pentatope (p0-p1-p2-p3-p4) of the topology
      result = GEO::PCK::side5_SOS(zi,site_[0],site_[1],site_[2],site_[3],zj,simplex_[0],simplex_[1],simplex_[2],simplex_[3],simplex_[4],dim_);
      avro_assert( result!=GEO::ZERO );
      break;
    default:
      printf("more predicates needed for %d-simplices!",(int)simplex_.size());
      avro_implement;
  }

  avro_assert( result!=GEO::ZERO );

  return result;
}

void
Vertex::print( const std::string& pre , const bool symbolic ) const
{
	printf("%s: (",pre.c_str());
	for (coord_t d=0;d<dim_;d++)
		printf(" %3.2e",x_[d]);
	printf(" )");
	if (symbolic)
	{
		// print the bisectors
		printf(" - [");
		for (index_t k=0;k<bisector_.size();k++)
			printf(" b%d ",bisector_[k]);
		printf("]");

		// print the simplex vertex information
		printf(" - [");
		for (index_t k=0;k<simplex_.size();k++)
			printf(" s%p ",(void*)simplex_[k]);
		printf("]");
		// print the topologies
		// TODO

		// print the sites
		printf(" - [");
		for (index_t k=0;k<site_.size();k++)
			printf(" z%p ",(void*)site_[k]);
		printf("]");
	}
	printf("\n");
}

} // delaunay

} // avro
