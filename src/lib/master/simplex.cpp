#include "common/tools.h"

#include "master/quadrature.h"
#include "master/simplex.h"

#include "mesh/topology.h"

namespace luna
{

Simplex::Simplex( const Topology<Simplex>& topology , const coord_t order ) :
  Simplex(topology.number(),order)
{}

void
Simplex::precalculate()
{
  // save the vertices
  // TODO

  // save the edges
  for (index_t k=0;k<number_+1;k++)
  for (index_t i=k+1;i<number_+1;i++)
  {
    this->edges_.push_back(k);
    this->edges_.push_back(i);
  }

  // triangulate the reference element
  // TODO

  // precalculate the shape function values at the quadrature points

}

void
Simplex::get_facet_vertices( const index_t* v , index_t nv , index_t ifacet , Element& f ) const
{
  f.indices.resize(f.dim+1);

  if (f.dim==0)
  {
    f.indices[0] = get_vertex(v,nv,ifacet);
    return;
  }
  if (f.dim==number_)
  {
    luna_assert( ifacet==0 );
    for (index_t j=0;j<nv;j++)
      f.indices[j] = v[j];
    if (f.sorted)
      std::sort( f.indices.begin() , f.indices.end() );
    return;
  }

  if (f.dim==1)
    get_edge( v,nv,ifacet,f.indices.data() );
  else if (f.dim==2)
    get_triangle( v,nv,ifacet,f.indices.data() );
  else if (f.dim==number_-1)
    get_facet_vertices( v,nv,ifacet,f.indices );
  else
    luna_implement;

  if (f.sorted)
    std::sort( f.indices.begin() , f.indices.end() );
}

index_t
Simplex::get_index( index_t dim , index_t ifacet , index_t ilocal ) const
{
  if (dim==0) return ifacet;
  if (dim==1) return number_+1 + ifacet*nb_interior(dim) + ilocal;
  if (dim==2) return number_+1 + nb_edges()*nb_interior(1) + nb_interior(2)*ifacet + ilocal;
  if (dim==3) return number_+1 + nb_edges()*nb_interior(1) + nb_interior(2)*nb_facets(2) + nb_interior(3)*ifacet + ilocal;
  luna_assert_not_reached;
  return nb_basis()+1;
}

index_t
Simplex::get_vertex( const index_t* v , index_t nv , index_t ivertex ) const
{
  return v[ivertex]; // vertices always stored at the beginning
}

void
Simplex::get_edge( const index_t* v , index_t nv , index_t iedge , index_t* e ) const
{
  index_t p0;
  index_t p1;

  luna_assert( number_>=1 );
  if (iedge==0)
  {
    p0 = get_vertex(v,nv,0);
    p1 = get_vertex(v,nv,1);
  }
  if (number_==1) goto sort;

  if (iedge==1)
  {
    p0 = get_vertex(v,nv,1);
    p1 = get_vertex(v,nv,2);
  }
  if (iedge==2)
  {
    p0 = get_vertex(v,nv,0);
    p1 = get_vertex(v,nv,2);
  }
  if (number_==2) goto sort;

  if (iedge==3)
  {
    p0 = get_vertex(v,nv,0);
    p1 = get_vertex(v,nv,3);
  }
  if (iedge==4)
  {
    p0 = get_vertex(v,nv,1);
    p1 = get_vertex(v,nv,3);
  }
  if (iedge==5)
  {
    p0 = get_vertex(v,nv,2);
    p1 = get_vertex(v,nv,3);
  }
  if (number_==3) goto sort;

  if (number_==4) luna_implement;

sort:
  if (p0<p1)
  {
    e[0] = p0;
    e[1] = p1;
  }
  else
  {
    e[0] = p1;
    e[1] = p0;
  }
}

void
Simplex::get_triangle( const index_t* v , index_t nv , index_t itriangle, index_t* t ) const
{
  // this is incorrect
  luna_implement;
  index_t p0 = get_vertex( v , nv ,  itriangle   %3 );
  index_t p1 = get_vertex( v , nv , (itriangle+1)%3 );
  index_t p2 = get_vertex( v , nv , (itriangle+2)%3 );
  t[0] = p0;
  t[1] = p1;
  t[2] = p2;
}

void
Simplex::get_facet_vertices( const index_t* v , index_t nv , index_t ifacet, std::vector<index_t>& f ) const
{
  f.resize( number_+1 );
  for (coord_t j=0;j<number_+1;j++)
    f[j] = v[j];
  f.erase( f.begin()+ifacet );
}

void
Simplex::get_edges( const index_t* v , const index_t nv , std::vector<index_t>& ek ) const
{
  // retrieve all edges of the simplex
  ek.resize( edges_.size() );
  for (index_t j=0;j<edges_.size();j++)
    ek[j] = v[edges_[j]];
}

} // luna
