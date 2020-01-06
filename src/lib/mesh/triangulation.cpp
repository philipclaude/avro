#include "master/master.h"
#include "master/polytope.h"
#include "master/simplex.h"

#include "mesh/points.h"
#include "mesh/triangulation.h"

#include "numerics/geometry.h"

#include <algorithm>
#include <set>

namespace avro
{

template<typename type>
Triangulation<type>::Triangulation( const Topology<type>& topology ) :
  TriangulationBase(topology.number(),topology.points().dim()),
  topology_(topology)
{
  // create topologies to store the lower-dimensional simplices
  for (index_t k=0;k<=number_;k++)
  {
    add_child( std::make_shared<Topology<Simplex>>(points_,k) );
    elements_.push_back( std::map<Element,index_t>() );
  }
}

template<typename type>
index_t
Triangulation<type>::add_simplex( index_t number , const index_t* v )
{
  std::vector<index_t> simplex(v,v+number+1);
  std::sort( simplex.begin() , simplex.end() );
  Element element;
  element.dim     = number;
  element.indices = simplex;
  element.sorted  = true;
  if (elements_[number].find(element)==elements_[number].end())
  {
    elements_[number].insert( {element,child_[number]->nb()} );
    child_[number]->add( v , number+1 );
  }
  return elements_[number].at(element);
}

template<typename type>
index_t
Triangulation<type>::add_point( coord_t number , const index_t* v , index_t nv )
{
  std::vector<index_t> polytope(v,v+nv);
  std::sort( polytope.begin() , polytope.end() );
  Element element;
  element.dim     = number;
  element.indices = polytope;
  element.sorted  = true;
  if (centroids_.find(element)==centroids_.end())
  {
    centroids_.insert( {element,points_.nb()} );
    centroid2dim_.insert( {points_.nb(),number} );
    std::vector<real_t> xc( points_.dim() );
    numerics::centroid( v , nv , topology_.points() , xc );
    points_.create( xc.data() );
  }
  else printf("retrieving point!!\n");
  return centroids_.at(element);
}

template<typename type>
void
Triangulation<type>::get_triangles( std::vector<index_t>& triangles ) const
{
  for (index_t k=0;k<child_[2]->nb();k++)
  {
    bool skip = false;
    for (index_t j=0;j<3;j++)
    {
      index_t idx = (*child_[2])(k,j);
      if (centroid2dim_.find(idx)==centroid2dim_.end()) continue;
      if (centroid2dim_.at(idx)!=2)
      {
        skip = true;
        break;
      }
    }
    if (skip) continue;
    for (index_t j=0;j<3;j++)
      triangles.push_back( (*child_[2])(k,j) );
  }
}

template<>
void
Triangulation<Simplex>::extract()
{
  std::vector<index_t> tk;
  std::set<std::string> MAP;
  std::vector<index_t> triangle(3);
  std::string s;

  for (index_t k=0;k<topology_.points().nb();k++)
    points_.create( topology_.points()[k] );

  // loop through all the cells
  for (index_t k=0;k<topology_.nb();k++)
  {
    if (topology_.ghost(k)) continue;

    // get the edges of this cell
    topology_.master().get_triangles( topology_(k) , topology_.nv(k) , tk );

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
        add_simplex( 2 , triangle.data() );
      }
    }
  }
}

template<>
void
Triangulation<Polytope>::extract()
{
  for (index_t k=0;k<topology_.points().nb();k++)
    points_.create( topology_.points()[k] );

  // loop through the cells
  for (index_t k=0;k<topology_.nb();k++)
  {
    if (topology_.ghost(k)) continue;

    // ask the master to triangulate, points will be added to points stored
    // in the triangulation object upon triangulation by the master
    topology_.master().triangulate( topology_(k) , topology_.nv(k) , *this );
  }
}

template class Triangulation<Simplex>;
template class Triangulation<Polytope>;

} // avro
