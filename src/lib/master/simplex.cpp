#include "common/tools.h"

#include "master/quadrature.h"
#include "master/simplex.h"

#include "mesh/points.h"
#include "mesh/topology.h"

#include "numerics/geometry.h"
#include "numerics/matrix.h"

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

real_t
Simplex::closest( const Points& x , const index_t* v , const index_t nv , const real_t* p , std::vector<real_t>& y ) const
{
  // for the point p, compute the closest barycentric point via least-squares
  // i.e. minimize || x(alpha)-p ||^2 s.t. sum(alpha) = 1 and all(alpha) > 0
  // the minimizer satisfiers
  // [V'*V, 1 ; 1' , 0] * [alpha;lambda] = [V'*p;1]
  // where V is the matrix of vertex coordinates and lambda is the
  // lagrange multiplier
  //
  //avro_assert_msg( nv==index_t(number+1) ,  );
  //avro_assert( alpha.size()==index_t(number+1) );

  // set the matrix of vertex coordinates
  #if 1
  luna_implement
  #else
  numerics::densMat<real> V(x.dim(),nv);
  for (index_t k=0;k<nv;k++)
  for (index_t j=0;j<x.dim();j++)
    V(j,k) = x[ v[k] ][j];
  const numerics::densMat<real> Vt = V.transpose();

  // set the least-squares portion of the system
  numerics::densMat<real> A = Vt*V;
  numerics::densMat<real> B(nv+1,nv+1);
  for (index_t i=0;i<nv;i++)
  for (index_t j=0;j<nv;j++)
    B(i,j) = A(i,j);

  // include the constraint in the full system
  for (index_t i=0;i<nv;i++)
  {
    B(i,nv) = 1;
    B(nv,i) = 1;
  }

  // set the right-hand side
  std::vector<real> b(nv+1,0.);
  Vt.multiply( p , b.data() );
  b[nv] = 1.0;

  // solve the system
  B.solve(b.data());

  // set the barycentric coordinates ignoring the last entry
  // because it holds the lagrange multiplier
  std::vector<real> alpha(nv);
  for (index_t i=0;i<nv;i++)
    alpha[i] = b[i];

  // check if the barycentric coordinates are in the simplex
  bool outside = false;
  for (index_t k=0;k<nv;k++)
  {
    if (alpha[k]<0)
    {
      outside = true;
      break;
    }
  }
  if (!outside)
  {
    // compute the actual coordinates and return the squared-distance
    V.multiply(alpha.data(),y.data());
    return geometrics::distance2(y.data(),p,x.dim());
  }

  // the projected point is outside the simplex
  // iterate through the facets
  std::vector<index_t> s(v,v+nv);
  std::vector<index_t> f;
  for (index_t j=0;j<nv;j++)
  {
    f = s;
    f.erase(f.begin()+j);

    std::vector<real> yf(x.dim());
    real df = closestBarycentric(x,f.data(),f.size(),p,yf,distance2);
    if (df<distance2 || distance2<0)
    {
      // the point minimizes the distance so set the point and save
      // the new squared-distance
      y = yf;
      distance2 = df;
    }
  }

  // no point was found that is both inside the simplex and minimizes the input squared-distance
  return distance2;
  #endif
}

static real_t
simplex_volume( const std::vector<const real_t*>& x , const coord_t dim )
{
  if (x.size()==3 && dim==2) return numerics::volume2(x);
  if (x.size()==4 && dim==3) return numerics::volume3(x);
  if (x.size()==5 && dim==4) return numerics::volume4(x);
  return numerics::volume_nd(x,dim);
}

real_t
Simplex::volume( const Points& points , const index_t* v , index_t nv ) const
{
  luna_assert( index_t(number_+1)==nv );
  std::vector<const real_t*> x(nv);
  for (index_t k=0;k<nv;k++)
    x[k] = points[v[k]];
  return simplex_volume(x,points.dim());
}

index_t
Simplex::edge( index_t k , index_t i ) const
{
  luna_assert( i <= 2 );
  return edges_[2*k+i];
}

} // luna
