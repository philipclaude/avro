//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//

#include "element/element.h"
#include "element/polytope.h"
#include "element/simplex.h"

#include "mesh/points.h"
#include "mesh/decomposition.h"

#include "numerics/geometry.h"

#include <algorithm>
#include <set>

namespace avro
{

template<typename type>
SimplicialDecomposition<type>::SimplicialDecomposition( const Topology<type>& topology ) :
  SimplicialDecompositionBase(topology.number(),topology.points().dim()),
  topology_(topology)
{
  // create topologies to store the lower-dimensional simplices
  for (index_t k=0;k<=number_;k++)
  {
    add_child( std::make_shared<Topology<Simplex>>(points_,k) );
    elements_.push_back( std::map<ElementIndices,index_t>() );
    parents_.push_back( std::map<index_t,std::vector<index_t>>() );
  }

  // copy all the points and set them as native
  for (index_t k=0;k<topology_.points().nb();k++)
  {
    points_.create( topology_.points()[k] );
    native_.push_back(true);
  }

  reference_coordinates_.set_layout(TableLayout_Jagged);
  point_parents_.resize( topology_.points().nb() );
  std::vector<bool> visited( topology_.points().nb() , false );
  for (index_t k=0;k<topology_.nb();k++)
  {
    for (index_t j=0;j<topology_.nv(k);j++)
    {
      if (visited[topology_(k,j)]) continue;
      visited[topology(k,j)] = true;

      std::vector<real_t> alpha(topology_.nv(k),0.0);
      alpha[j] = 1.0;

      avro_assert( topology_(k,j) < topology_.points().nb() );

      reference_coordinates_.add( alpha.data() , alpha.size() );
      point_parents_[ topology_(k,j) ] = k;
    }
  }
}

template<typename type>
index_t
SimplicialDecomposition<type>::add_simplex( index_t number , const index_t* v , index_t parent )
{
  std::vector<index_t> simplex(v,v+number+1);
  std::sort( simplex.begin() , simplex.end() );
  ElementIndices element;
  element.dim     = number;
  element.indices = simplex;
  element.sorted  = true;
  std::map<ElementIndices,index_t>::iterator it = elements_[number].find(element);
  if (it==elements_[number].end())
  {
    elements_[number].insert( { element,child_[number]->nb() } );
    parents_[number].insert( {child_[number]->nb(),{}} );
    child_[number]->add( v , number+1 );
  }
  index_t idx = elements_[number].at(element);
  parents_[number][idx].push_back(parent);
  return elements_[number].at(element);
}

template<typename type>
index_t
SimplicialDecomposition<type>::add_point( coord_t number , const index_t* v , index_t nv , index_t parent )
{
  std::vector<index_t> polytope(v,v+nv);
  std::sort( polytope.begin() , polytope.end() );
  ElementIndices element;
  element.dim     = number;
  element.indices = polytope;
  element.sorted  = true;
  if (true) //centroids_.find(element)==centroids_.end()) // TODO this caused tests to fail (with map.at) on wazowski but not sure why...
  {
    centroids_.insert( {element,points_.nb()} );
    centroid2dim_.insert( {points_.nb(),number} );
    std::vector<real_t> xc( points_.dim() );
    numerics::centroid( v , nv , topology_.points() , xc );
    points_.create( xc.data() );
    native_.push_back(false);
    point_parents_.push_back(parent);

    std::vector<real_t> alpha( topology_.nv(parent) , 0 );
    for (index_t j=0;j<topology_.nv(parent);j++)
    {
      for (index_t i=0;i<nv;i++)
      {
        if (v[i]==topology_(parent,j))
        {
          alpha[j] = 1./nv;
          break;
        }
      }
    }
    reference_coordinates_.add( alpha.data() , alpha.size() );
  }
  return points_.nb()-1;//centroids_.at(element);
}

template<typename type>
void
SimplicialDecomposition<type>::get_simplices( coord_t number , std::vector<index_t>& simplices , std::vector<index_t>& parents ) const
{
  const Topology<Simplex>& s = *child_[number];

  for (index_t k=0;k<s.nb();k++)
  {
    if (topology_.number()==2)
      avro_assert_msg( parents_[number].at(k).size()==1 , "nb_parents of elem %lu = %lu" , k , parents_[number].at(k).size() ); // all triangles have cardinality 1
    //if (topology_.number()==3)
    //  avro_assert_msg( parents_[number].at(k).size()<=2 , "cardinality = %lu" , parents_[number].at(k).size() ); // boundary facets have cardinality 1, interior 2

    //if (topology_.number()==3 && parents_[number].at(k).size()!=1) continue;

    for (coord_t j=0;j<number+1;j++)
      simplices.push_back(s(k,j));

    parents.push_back( parents_[number].at(k)[0] );
  }
}

template<>
void
SimplicialDecomposition<Simplex>::extract()
{
  std::vector<index_t> tk;
  std::set<std::string> MAP;
  std::vector<index_t> triangle(3);
  std::string s;

  // loop through all the cells
  for (index_t k=0;k<topology_.nb();k++)
  {
    if (topology_.ghost(k)) continue;

    // get the edges of this cell
    topology_.element().get_triangles( topology_(k) , topology_.nv(k) , tk );

    // add the edges
    for (index_t j=0;j<tk.size()/3;j++)
    {
      index_t p0 = tk[3*j];
      index_t p1 = tk[3*j+1];
      index_t p2 = tk[3*j+2];

      if (p0<topology_.points().nb_ghost() ||
          p1<topology_.points().nb_ghost() ||
          p2<topology_.points().nb_ghost())
        continue;

      triangle[0] = p0;
      triangle[1] = p1;
      triangle[2] = p2;
      s = unique_label(triangle);

      if (MAP.find(s)==MAP.end())
      {
        MAP.insert(s);
        add_simplex( 2 , triangle.data() , k );
      }
    }
  }
  if (number_==topology_.points().dim())
    child(number_).orient();
}

template<>
void
SimplicialDecomposition<Polytope>::extract()
{
  // loop through the cells
  for (index_t k=0;k<topology_.nb();k++)
  {
    if (topology_.ghost(k)) continue;

    // ask the  to triangulate, points will be added to points stored
    // in the SimplicialDecomposition object upon decomposition by the
    std::set<int> h;
    topology_.element().triangulate( topology_(k) , topology_.nv(k) , *this , k , h );
    avro_assert( h.size()==0 );
  }
  if (number_==topology_.points().dim())
    child(number_).orient();
}

template class SimplicialDecomposition<Simplex>;
template class SimplicialDecomposition<Polytope>;

} // avro
