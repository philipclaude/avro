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
    parents_.push_back( std::map<index_t,std::vector<index_t>>() );
  }

  // copy all the points and set them as native
  for (index_t k=0;k<topology_.points().nb();k++)
  {
    points_.create( topology_.points()[k] );
    native_.push_back(true);
  }

  barycentric_.set_layout(TableLayout_Jagged);

  point2element_.resize( topology_.points().nb() );
  std::vector<bool> visited( topology_.points().nb() , false );
  for (index_t k=0;k<topology_.nb();k++)
  {
    for (index_t j=0;j<topology_.nv(k);j++)
    {
      if (visited[topology_(k,j)]) continue;
      visited[topology(k,j)] = true;

      Element element;
      element.indices.assign( topology_(k) , topology_(k)+topology_.nv(k) );
      element.dim = topology_.number();
      element.sorted = false;
      point2element_[topology(k,j)] = element;
      std::vector<real_t> alpha(topology_.nv(k),0.0);
      alpha[j] = 1.0;
      barycentric_.add( alpha.data() , alpha.size() );
    }
  }
}

template<typename type>
index_t
Triangulation<type>::add_simplex( index_t number , const index_t* v , index_t parent )
{
  std::vector<index_t> simplex(v,v+number+1);
  std::sort( simplex.begin() , simplex.end() );
  Element element;
  element.dim     = number;
  element.indices = simplex;
  element.sorted  = true;
  std::map<Element,index_t>::iterator it = elements_[number].find(element);
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
    native_.push_back(false);
    point2element_.push_back( element );

    std::vector<real_t> alpha( nv , 1./nv );
    barycentric_.add( alpha.data() , alpha.size() );
  }
  return centroids_.at(element);
}

template<typename type>
void
Triangulation<type>::get_simplices( coord_t number , std::vector<index_t>& simplices , std::vector<index_t>& parents ) const
{
  const Topology<Simplex>& s = *child_[number];

  for (index_t k=0;k<s.nb();k++)
  {
    bool skip = false;
    for (index_t j=0;j<number+1;j++)
    {
      index_t idx = s(k,j);
      if (centroid2dim_.find(idx)==centroid2dim_.end()) continue;
      if (centroid2dim_.at(idx)!=number)
      {
        skip = true;
        break;
      }
    }
    if (skip) continue;
    for (index_t j=0;j<number+1;j++)
      simplices.push_back(s(k,j));

    if (topology_.number()==2)
      avro_assert( parents_[number].at(k).size()==1 ); // all triangles have cardinality 1
    if (topology_.number()==3)
      avro_assert( parents_[number].at(k).size()<=2 ); // boundary facets have cardinality 1, interior 2

    parents.push_back( parents_[number].at(k)[0] );
  }
}

template<>
void
Triangulation<Polytope>::get_barycentric( coord_t number , Table<index_t>& dof , Table<real_t>& barycentric ) const
{
  dof.set_layout( TableLayout_Jagged );
  barycentric.set_layout( TableLayout_Jagged );

  avro_assert( point2element_.size() == points_.nb() );
  avro_assert( barycentric_.nb() == points_.nb() );

  // create interpolation data for each point
  for (index_t k=0;k<points_.nb();k++)
  {
    // add the dof data for this point
    const Element& element = point2element_[k];
    print_inline(element.indices);
    avro_assert( element.dim == number );

    dof.add( element.indices.data() , element.indices.size() );
    barycentric.add( barycentric_(k) , barycentric_.nv(k) );

    std::vector<real_t> alpha( barycentric_(k) , barycentric_(k)+barycentric_.nv(k) );
    print_inline(alpha);
  }
}

template<>
void
Triangulation<Simplex>::get_barycentric( coord_t number , Table<index_t>& dof , Table<real_t>& barycentric ) const
{
  avro_implement;
}

template<>
void
Triangulation<Simplex>::extract()
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
        add_simplex( 2 , triangle.data() , k );
      }
    }
  }
  child(number_).orient();
}

template<>
void
Triangulation<Polytope>::extract()
{
  // loop through the cells
  for (index_t k=0;k<topology_.nb();k++)
  {
    if (topology_.ghost(k)) continue;

    // ask the master to triangulate, points will be added to points stored
    // in the triangulation object upon triangulation by the master
    topology_.master().triangulate( topology_(k) , topology_.nv(k) , *this , k );
  }
  child(number_).orient();
}

template class Triangulation<Simplex>;
template class Triangulation<Polytope>;

} // avro
