#include "master/polytope.h"
#include "master/simplex.h"

#include "mesh/points.h"
#include "mesh/triangulation.h"

#include <algorithm>
#include <set>

namespace avro
{

template<typename type>
Triangulation<type>::Triangulation( const Topology<type>& topology ) :
  TriangulationBase(points_,topology.number()),
  topology_(topology),
  points_(topology.points().dim())
{}

template<>
void
Triangulation<Simplex>::extract()
{
  std::vector<index_t> tk;
  std::set<std::string> MAP;
  std::vector<index_t> triangle(3);
  std::string s;

  // loop through all the cells
  for (index_t k=0;k<nb();k++)
  {
    if (ghost(k)) continue;

    const index_t* v0 = topology_(k);

    // get the edges of this cell
    topology_.master().get_triangles( v0 , nv(k) , tk );

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
        add( triangle.data() , triangle.size() );
      }
    }
  }
}

template<>
void
Triangulation<Polytope>::extract()
{
  if (number_==0)
  {
    // add the simplices
    for (index_t k=0;k<topology_.nb();k++)
      add( topology_(k) , 1 );
  }

  // loop through the cells
  for (index_t k=0;k<topology_.nb();k++)
  {

    if (topology_.ghost(k)) continue;

    // record the number of vertices and simplices we currently have
    index_t nv0 = points_.nb();
    index_t nt0 = this->nb();

    // ask the master to triangulate -- it needs these vertices to compute the
    // geometry of the vertices it creates
    topology_.master().triangulate( topology_(k) , topology_.nv(k) , *this );

    index_t nt = this->nb();

    // loop through the created simplices
    for (index_t j=nt0;j<nt;j++)
    {
      avro_assert( this->nv(j)==index_t(this->number()+1) );
      for (index_t i=0;i<this->nv(j);i++)
      {
        // add the parent data for this vertex
        parents_.push_back( k );
      }
    }

    // add these vertices to the vertex list for the interpolator (graphics)
    for (index_t j=nv0;j<topology_.points().nb();j++)
      points_.create( topology_.points()[k] );
  }
}

template class Triangulation<Simplex>;
template class Triangulation<Polytope>;

} // avro
