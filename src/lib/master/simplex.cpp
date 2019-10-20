#include "common/tools.h"

#include "master/quadrature.h"
#include "master/simplex.h"

#include "mesh/topology.h"

namespace ursa
{

template<typename Basis>
SimplexBase<Basis>::SimplexBase( const Topology<Simplex<Basis>>& topology , const coord_t order ) :
  SimplexBase(topology.number(),order)
{}

template<typename Basis>
void
SimplexBase<Basis>::get_facet_vertices( const index_t* v , index_t nv , index_t ifacet , Element& f ) const
{
  f.indices.resize(f.dim+1);

  if (f.dim==0)
  {
    f.indices[0] = get_vertex(v,nv,ifacet);
    return;
  }
  if (f.dim==number_)
  {
    ursa_assert( ifacet==0 );
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
    ursa_implement;

  if (f.sorted)
    std::sort( f.indices.begin() , f.indices.end() );
}

template<typename Basis>
index_t
SimplexBase<Basis>::get_index( index_t dim , index_t ifacet , index_t ilocal ) const
{
  if (dim==0) return ifacet;
  if (dim==1) return number_+1 + ifacet*nb_interior(dim) + ilocal;
  if (dim==2) return number_+1 + nb_edges()*nb_interior(1) + nb_interior(2)*ifacet + ilocal;
  if (dim==3) return number_+1 + nb_edges()*nb_interior(1) + nb_interior(2)*nb_facets(2) + nb_interior(3)*ifacet + ilocal;
  ursa_assert_not_reached;
  return nb_basis()+1;
}

template<typename Basis>
index_t
SimplexBase<Basis>::get_vertex( const index_t* v , index_t nv , index_t ivertex ) const
{
  // maybe this only specializes for lagrange?
  return v[ivertex]; // vertices always stored at the beginning
}

template<typename Basis>
void
SimplexBase<Basis>::get_edge( const index_t* v , index_t nv , index_t iedge , index_t* e ) const
{
  index_t p0;
  index_t p1;

  ursa_assert( number_>=1 );
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

  if (number_==4) ursa_implement;

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

template<typename Basis>
void
SimplexBase<Basis>::get_edges( const index_t* v , const index_t nv , std::vector<index_t>& ek ) const
{
  ek.resize( edges_.size() );
  for (index_t j=0;j<edges_.size();j++)
    ek[j] = v[edges_[j]];
}

template<typename Basis>
void
SimplexBase<Basis>::get_triangle( const index_t* v , index_t nv , index_t itriangle, index_t* t ) const
{
  index_t p0 = get_vertex( v , nv ,  itriangle   %3 );
  index_t p1 = get_vertex( v , nv , (itriangle+1)%3 );
  index_t p2 = get_vertex( v , nv , (itriangle+2)%3 );
  t[0] = p0;
  t[1] = p1;
  t[2] = p2;
}

template<typename Basis>
void
SimplexBase<Basis>::get_facet_vertices( const index_t* v , index_t nv , index_t ifacet, std::vector<index_t>& f ) const
{
  f.resize( number_+1 );
  for (coord_t j=0;j<number_+1;j++)
    f[j] = v[j];
  f.erase( f.begin()+ifacet );
}

template<typename Basis>
void
SimplexBase<Basis>::loadQuadrature( Quadrature& quadrature )
{
  quadrature.retrieve(xquad_,wquad_);
}

template class SimplexBase<Lagrange>;
template class SimplexBase<Bezier>;


} // ursa
