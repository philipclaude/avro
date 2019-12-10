#include "mesh/points.h"
#include "mesh/triangulation.h"

namespace luna
{

template<typename type>
Triangulation<type>::Triangulation( const Topology<type>& topology ) :
  Topology<Simplex>(points_,topology.number(),1),
  topology_(topology)
{}

template<typename type>
void
Triangulation<type>::extract()
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
    topology_.master().triangulate( number() , *this , points_ ,
                                    topology_(k) , topology_.nv(k) );

    index_t nt = this->nb();

    // loop through the created simplices
    for (index_t j=nt0;j<nt;j++)
    {
      luna_assert( this->nv(j)==index_t(this->number()+1) );
      for (index_t i=0;i<this->nv(j);i++)
      {
        // add the parent data for this vertex
        parents_.push_back( k );
      }
    }

    // add these vertices to the vertex list for the interpolator (graphics)
    //for (index_t j=nv0;j<points_.nb();j++)
    //  points_.create( points_.get(j) );
  }
}

template class Triangulation<Simplex>;

} // luna
