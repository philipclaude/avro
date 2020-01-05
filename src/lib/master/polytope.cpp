#include "common/error.h"
#include "common/set.h"
#include "common/tools.h"

#include "master/polytope.h"

#include "mesh/points.h"
#include "mesh/topology.h"
#include "mesh/triangulation.h"

#include "numerics/geometry.h"

namespace avro
{

Polytope::Polytope( coord_t number , coord_t order , const Table<int>& incidence ) :
  Master(number,order),
  simplex_(number,order),
  incidence_(incidence)
{
  avro_assert_msg( order==1 , "not supported..." );
}

Polytope::Polytope( Topology<Polytope>& topology , const coord_t order , const Table<int>& incidence ) :
  Master(topology.number(),order),
  simplex_(topology.number(),order),
  incidence_(incidence)
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

void
Polytope::triangulate( const index_t* v , index_t nv , Triangulation<Polytope>& triangulation ) const
{
  const Topology<Polytope>& topology = triangulation.topology();
  const Points& points = topology.points();

  // construct a lower dimensional polytope with the same vertex facet matrix
  Polytope facetope(number_-1,order_,incidence_);

  if (triangulation.number()<=1) return;

  // compute the geometry of the vertex we will create
  std::vector<real_t> xc( points.dim() );
  numerics::centroid( v , nv , points , xc );

  // save the id of the vertex we create
  int id = points.nb();

  // create the new vertex with number 'id'
  triangulation.points().create( xc );

  // get the hrep of this cell
  std::vector<int> facets;
  hrep( v ,  nv , facets );

  // loop through each facet
  for (index_t j=0;j<facets.size();j++)
  {
    // get the points with this bisector
    std::vector<index_t> vf;
    vrep( v , nv , facets[j] , vf );

    // this also defines a polytope which we need to triangulate!
    Table<index_t> tf(TableLayout_Jagged);

    if (vf.size()==index_t(number_))
    {
      // add the vrep to the triangulation
      tf.add( vf.data() , vf.size() );
    }

    // triangulate the lower dimensional polytope
    else facetope.triangulate( number_-1 , tf , triangulation.points() , vf.data() , vf.size() );

    // now connect the resulting simplices to the created vertex
    if (number_==3 && triangulation.number()==2)
    {
      // this was a 3d polytope but we want the bounding triangles, likely to plot
      // do not connect the facets, simply store them
      for (index_t k=0;k<tf.nb();k++)
      {
        std::vector<index_t> tk(3);
        for (index_t i=0;i<tf.nv(k);i++)
          tk[i] = tf(k,i);
        triangulation.add(tk.data(),tk.size());
      }
    }
    else
    {
      for (index_t k=0;k<tf.nb();k++)
      {
        avro_assert( tf.nv(k)==index_t(number_) );
        std::vector<index_t> tk(number_+1);
        for (index_t i=0;i<tf.nv(k);i++)
          tk[i] = tf(k,i);
        tk[ number_ ] = id; // the added vertex

        // add the new simplex
        triangulation.add( tk.data() , tk.size() );
      }
    }
  }
}

void
Polytope::triangulate( coord_t number , Table<index_t>& simplices , Points& points , const index_t* v , index_t nv ) const
{
  // construct a lower dimensional polytope with the same vertex facet matrix
  Polytope facetope(number_-1,order_,incidence_);

  // nothing to triangulate if we have reached points or edges
  if (number<=1) return;

  // compute the coordinates of the point we will create
  std::vector<real_t> xc( points.dim() );
  numerics::centroid( v , nv , points , xc );

  // save the id of the point we create
  int id = points.nb();

  // create the new vertex with number 'id'
  points.create( xc );

  // get the hrep of this cell
  std::vector<int> facets;
  hrep( v ,  nv , facets );

  // loop through each facet
  for (index_t j=0;j<facets.size();j++)
  {
    // get the points with this bisector
    std::vector<index_t> vf;
    vrep( v , nv , facets[j] , vf );

    // this also defines a polytope which we need to triangulate!
    Table<index_t> facets(TableLayout_Jagged);

    if (vf.size()==index_t(number_))
    {
      // add the vrep to the triangulation
      facets.add( vf.data() , vf.size() );
    }

    // triangulate the lower dimensional polytope
    else facetope.triangulate( number-1 , facets , points , vf.data() , vf.size() );

    // now connect the resulting simplices to the created vertex
    if (number_==3 && number==2)
    {
      // this was a 3d polytope but we want the bounding triangles, likely to plot
      // do not connect the facets, simply store them
      for (index_t k=0;k<facets.nb();k++)
      {
        std::vector<index_t> tk(3);
        for (index_t i=0;i<facets.nv(k);i++)
          tk[i] = facets(k,i);
        simplices.add(tk.data(),tk.size());
      }
    }
    else
    {
      for (index_t k=0;k<facets.nb();k++)
      {
        avro_assert( facets.nv(k)==index_t(number_) );
        std::vector<index_t> tk(number_+1);
        for (index_t i=0;i<facets.nv(k);i++)
          tk[i] = facets(k,i);
        tk[ number_ ] = id; // the added vertex

        // add the new simplex
        simplices.add( tk.data() , tk.size() );
      }
    }
  }
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
