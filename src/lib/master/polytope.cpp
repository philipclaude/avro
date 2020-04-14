#include "common/error.h"
#include "common/set.h"
#include "common/tools.h"

#include "master/polytope.h"

#include "mesh/decomposition.h"
#include "mesh/points.h"
#include "mesh/topology.h"

#include "numerics/geometry.h"

namespace avro
{

Polytope::Polytope( coord_t number , coord_t order , const Table<int>& incidence ) :
  Master(number,order),
  simplex_(number,order),
  incidence_(incidence),
  fullmesh_(false)
{
  avro_assert_msg( order==1 , "not supported..." );
}

Polytope::Polytope( Topology<Polytope>& topology , const coord_t order , const Table<int>& incidence ) :
  Master(topology.number(),order),
  simplex_(topology.number(),order),
  incidence_(incidence),
  fullmesh_(false)
{}


void
Polytope::get_edges( const index_t* v , index_t nv , std::vector<index_t>& ek ) const
{
  // given the vrep and the vertex-facet relations stored in facets_, construct the polytope edges
  ek.clear();

  for (index_t v0=0;v0<nv;v0++)
  {
    for (index_t v1=v0+1;v1<nv;v1++)
    {
      if ( is_edge(*(v+v0),*(v+v1)) )
      {
        ek.push_back( *(v+v0) );
        ek.push_back( *(v+v1) );
      }
    }
  }
}

real_t
Polytope::volume( const Points& points , const index_t* v , index_t nv ) const
{
  avro_implement;
  return -1.0;
}

void
Polytope::edges( const index_t* v , const index_t nv, std::vector<index_t>& e ) const
{
  // given the vrep and the vertex-facet relations stored in facets_, construct the polytope edges
  e.clear();

  for (index_t v0=0;v0<nv;v0++)
  {
    for (index_t v1=v0+1;v1<nv;v1++)
    {
      if ( is_edge(*(v+v0),*(v+v1)) )
      {
        e.push_back( *(v+v0) );
        e.push_back( *(v+v1) );
      }
    }
  }
}

void
Polytope::hrep( const index_t* v , index_t nv , std::vector<int>& facets ) const
{
  facets.clear();

  // loop through the points of this polytope
  for (index_t k=0;k<nv;k++)
  {
    // accumulate the facets of this vertex
    index_t v0 = *(v +k);
    for (index_t j=0;j<incidence_.nv(v0);j++)
      facets.push_back( incidence_(v0,j) );
  }
  uniquify(facets);
}

void
Polytope::vrep( const index_t* v , index_t nv , const int facet , std::vector<index_t>& points ) const
{
  points.clear();

  // loop through the points of this polytope
  for (index_t k=0;k<nv;k++)
  {
    // check if the vertex-facet relations of this vertex contain the facet
    if (incidence_.has( *(v+k) , facet ))
      points.push_back( *(v+k) );
  }
}

std::vector<index_t>
Polytope::triangulate( const index_t* v , index_t nv , SimplicialDecomposition<Polytope>& decomposition , index_t parent , std::set<int>& h ) const
{
  std::vector<index_t> simplex_idx;
  if (number_<=1)
  {
    avro_assert_msg( nv == number_+1 , "nv = %lu , number = %u" , nv , number_ );
    index_t idx = decomposition.add_simplex( number_ , v , parent );
    simplex_idx.push_back(idx);
    return simplex_idx;
  }

  // construct a lower dimensional polytope with the same vertex facet matrix
  Polytope facetope(number_-1,order_,incidence_); // TODO save this in class to avoid constructing simplex every time

  // the triangulation keeps track of which points have been added
  // and from which facet points symbolically created that point
  index_t id = 0;
  if (nv>index_t(number_+1))
    id = decomposition.add_point( number_ , v , nv , parent );

  // get the hrep of this cell
  std::vector<int> facets;
  hrep( v ,  nv , facets );

  // loop through each facet
  std::vector<index_t> facet_idx;
  for (index_t j=0;j<facets.size();j++)
  {
    // skip redundant bisectors
    if (h.find(facets[j])!=h.end()) continue;

    // get the points with this bisector
    std::vector<index_t> vf;
    vrep( v , nv , facets[j] , vf );

    // triangulate the lower dimensional polytope
    h.insert( facets[j] );
    std::vector<index_t> idx = facetope.triangulate( vf.data() , vf.size() , decomposition , parent , h );
    h.erase( facets[j] );

    // add the lower-dimensional facet identifiers
    for (index_t k=0;k<idx.size();k++)
      facet_idx.push_back(idx[k]);
  }

  // loop through all the simplicial facets
  if (nv==index_t(number_+1))
  {
    simplex_idx.push_back( decomposition.add_simplex( number_ , v , parent ) );
  }
  else
  {
    avro_assert( id > 0 );
    avro_assert( nv > number_ );
    for (index_t k=0;k<facet_idx.size();k++)
    {
      std::vector<index_t> simplex(number_+1);
      for (index_t i=0;i<number_;i++)
        simplex[i] = decomposition.child(number_-1)(facet_idx[k],i);
      simplex[number_] = id; // the added vertex

      // add the new simplex
      index_t idx = decomposition.add_simplex( number_ , simplex.data() , parent );
      simplex_idx.push_back(idx);
    }
  }

  return simplex_idx;
}

bool
Polytope::is_edge( index_t v0 , index_t v1 ) const
{
  // do the set intersection of the facets of each vertex
  // note the intersection assumes the sets to be sorted
  std::vector<int> facets;
  std::vector<int> f0,f1;
  f0 = incidence_.get(v0);
  f1 = incidence_.get(v1);
  Set::intersection( incidence_.get(v0) , incidence_.get(v1) , facets );
  if (fullmesh_)
    return facets.size()>=index_t(number_-1);
  return facets.size()==index_t(number_-1);
}

bool
Polytope::is_edge( const std::vector<int>& b0 , const std::vector<int>& b1 ) const
{
  // do the set intersection of the facets of each vertex
  // note the intersection assumes the sets to be sorted
  index_t count = 0;
  for (index_t k=0;k<b0.size();k++)
  {
    if (Set::contains(b0[k],b1)>-1)
      count++;
  }
  if (fullmesh_)
    return count>=index_t(number_-1);
  return count==index_t(number_-1);
}


} // avro
