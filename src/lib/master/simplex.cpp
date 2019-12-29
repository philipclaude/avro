#include "common/tools.h"

#include "master/quadrature.h"
#include "master/simplex.h"

#include "mesh/points.h"
#include "mesh/topology.h"

#include "numerics/geometry.h"
#include "numerics/linear_algebra.h"
#include "numerics/matrix.h"

namespace luma
{

Simplex::Simplex( const Topology<Simplex>& topology , const coord_t order ) :
  Simplex(topology.number(),order)
{}

Simplex::Simplex( const coord_t number , const coord_t order ) :
  Master(number,order,"simplex")
{
  precalculate();
}

index_t
Simplex::triangle( const index_t itriangle , const index_t inode ) const
{
  return triangles_[3*itriangle+inode];
}

std::vector<index_t>
Simplex::facet( const index_t j , const index_t i ) const
{
  std::vector<index_t> f(j+1);
  if (j==0) f[0] = i;
  else if (j==1)
  {
    f[0] = edge(i,0);
    f[1] = edge(i,1);
  }
  else if (j==2)
  {
    f[0] = triangle(i,0);
    f[1] = triangle(i,1);
    f[2] = triangle(i,2);
  }
  else if (j==index_t(number_-1))
  {
    index_t count = 0;
    for (index_t k=0;k<index_t(number_+1);k++)
    {
      if (k==i) continue;
      f[count++] = k;
    }
  }
  else
    luma_assert_not_reached;
  return f;
}

void
Simplex::precalculate()
{
  // save the vertices
  // TODO

  // save the edges
  for (index_t k=0;k<number_+1;k++)
  for (index_t i=k+1;i<number_+1;i++)
  {
    edges_.push_back(k);
    edges_.push_back(i);
  }

  // save the triangles
  for (coord_t i=0;i<number_+1;i++)
  {
    for (coord_t j=i+1;j<number_+1;j++)
    {
      for (coord_t k=j+1;k<number_+1;k++)
      {
        triangles_.push_back( i );
        triangles_.push_back( j );
        triangles_.push_back( k );
      }
    }
  }

  // triangulate the reference element
  // TODO

  // precalculate the shape function values at the quadrature points

}

void
Simplex::get_canonical_indices( const index_t* v , index_t nv , const Element& f , std::vector<index_t>& canonical ) const
{
  luma_assert_msg( order_==1 , "should this only be called for linear simplices?" );
  for (index_t k=0;k<f.indices.size();k++)
  {
    bool found = false;
    for (index_t j=0;j<nv;j++)
    {
      if (f.indices[k]==v[j])
      {
        canonical[k] = j;
        found = true;
        break;
      }
    }
    luma_assert( found );
  }
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
    luma_assert( ifacet==0 );
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
    luma_implement;

  if (f.sorted)
    std::sort( f.indices.begin() , f.indices.end() );
}

void
Simplex::facet( const index_t* v , const index_t j , std::vector<index_t>& f ) const
{
  index_t count = 0;
  for (index_t i=0;i<index_t(number_+1);i++)
  {
    if (i==j) continue;
    f[count++] = v[i];
  }
}

index_t
Simplex::get_vertex( const index_t* v , index_t nv , index_t ivertex ) const
{
  // not true for high-order simplices!!!
  return v[ivertex]; // vertices always stored at the beginning
}

void
Simplex::get_edge( const index_t* v , index_t nv , index_t iedge , index_t* e ) const
{
  index_t p0;
  index_t p1;

  luma_assert( number_>=1 );
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

  if (iedge==6)
  {
    p0 = get_vertex(v,nv,0);
    p1 = get_vertex(v,nv,4);
  }
  if (iedge==7)
  {
    p0 = get_vertex(v,nv,1);
    p1 = get_vertex(v,nv,4);
  }
  if (iedge==8)
  {
    p0 = get_vertex(v,nv,2);
    p1 = get_vertex(v,nv,4);
  }
  if (iedge==9)
  {
    p0 = get_vertex(v,nv,3);
    p1 = get_vertex(v,nv,4);
  }
  if (number_==4) goto sort;

  luma_implement;

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
  t[0] = v[triangles_[itriangle*3  ]];
  t[1] = v[triangles_[itriangle*3+1]];
  t[2] = v[triangles_[itriangle*3+2]];
  return;
  luma_implement;
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

void
Simplex::get_triangles( const index_t* v , const index_t nv , std::vector<index_t>& tk ) const
{
  // retrieve all edges of the simplex
  tk.resize( triangles_.size() );
  for (index_t j=0;j<triangles_.size();j++)
    tk[j] = v[triangles_[j]];
}

real_t
Simplex::closest( const Points& x , const index_t* v , const index_t nv , const real_t* p0 , std::vector<real_t>& y , real_t distance2 ) const
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
  #if 0
  luma_implement
  #else
  numerics::VectorD<real_t> p(x.dim(),p0);
  numerics::MatrixD<real_t> V(x.dim(),nv);
  for (index_t k=0;k<nv;k++)
  for (index_t j=0;j<x.dim();j++)
    V(j,k) = x[ v[k] ][j];
  //const numerics::densMat<real> Vt = V.transpose();

  // set the least-squares portion of the system
  //numerics::densMat<real> A = Vt*V;
  //numerics::densMat<real> B(nv+1,nv+1);
  numerics::SymMatrixD<real_t> A = numpack::Transpose(V)*V;
  numerics::MatrixD<real_t> B(nv+1,nv+1);
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
  //std::vector<real> b(nv+1,0.);
  //Vt.multiply( p , b.data() );
  numerics::VectorD<real_t> b0(nv+1);
  b0 = numpack::Transpose(V)*p;
  b0[nv] = 1.0;

  // solve the system
  numerics::VectorD<real_t> b(nv+1);
  b = numpack::DLA::InverseLU::Solve(B,b0);
  //B.solve(b.data());

  // set the barycentric coordinates ignoring the last entry
  // because it holds the lagrange multiplier
  //std::vector<real> alpha(nv);
  numerics::VectorD<real_t> alpha(nv);
  for (index_t i=0;i<nv;i++)
    alpha(i) = b[i];

  // check if the barycentric coordinates are in the simplex
  bool outside = false;
  for (index_t k=0;k<nv;k++)
  {
    if (alpha(k)<0)
    {
      outside = true;
      break;
    }
  }
  if (!outside)
  {
    // compute the actual coordinates and return the squared-distance
    numerics::VectorD<real_t> y0(x.dim());
    y0 = V*alpha;
    //V.multiply(alpha.data(),y.data());
    for (coord_t i=0;i<x.dim();i++)
      y[i] = y0(i);
    return numerics::distance2(y.data(),p0,x.dim());
  }

  // the projected point is outside the simplex
  // iterate through the facets
  std::vector<index_t> s(v,v+nv);
  std::vector<index_t> f;
  for (index_t j=0;j<nv;j++)
  {
    f = s;
    f.erase(f.begin()+j);

    std::vector<real_t> yf(x.dim());
    real_t df = closest(x,f.data(),f.size(),p0,yf,distance2);
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
  luma_assert_not_reached;
  return numerics::volume_nd(x,dim);
}

real_t
Simplex::volume( const Points& points , const index_t* v , index_t nv ) const
{
  luma_assert( index_t(number_+1)==nv );
  std::vector<const real_t*> x(nv);
  for (index_t k=0;k<nv;k++)
    x[k] = points[v[k]];
  return simplex_volume(x,points.dim());
}

index_t
Simplex::edge( index_t k , index_t i ) const
{
  luma_assert( i <= 2 );
  return edges_[2*k+i];
}

void
Simplex::triangulate( coord_t number , Topology<Simplex>& topology , Points& points , const index_t* v , index_t nv ) const
{
  luma_assert(nv==index_t(number+1));
  std::vector<index_t> tk(number+1);
  if (number==2)
  {
    // we were requested triangles
    for (index_t k=0;k<triangles_.size()/3;k++)
    {
      for (index_t j=0;j<3;j++)
        tk[j] = v[ triangles_[3*k+j] ];
      topology.add(tk.data(),tk.size());
    }
    return;
  }

  luma_assert(number==number_);
  for (index_t j=0;j<index_t(number+1);j++)
    tk[j] = v[j];
  topology.add(tk.data(),tk.size());
}

real_t
Simplex::jacobian( const std::vector<const real_t*>& x , coord_t dim ) const
{
  return simplex_volume(x,dim);
}

void
Simplex::jacobian( const std::vector<const real_t*>& xk , numerics::MatrixD<real_t>& J ) const
{
  luma_assert( J.m()==number_ );
  luma_assert( J.n()==number_ );
  luma_assert( order_==1 );
  for (index_t j=0;j<number_;j++)
  for (index_t d=0;d<number_;d++)
    J(d,j) = xk[j+1][d] -xk[0][d];
}

void
Simplex::jacobian( const index_t* v , index_t nv , const Points& points , numerics::MatrixD<real_t>& J ) const
{
  luma_assert( J.m()==number_ );
  luma_assert( J.n()==number_ );
  luma_assert( order_==1 );
  for (index_t j=0;j<number_;j++)
  for (index_t d=0;d<points.dim();d++)
    J(d,j) = points[v[j+1]][d] -points[v[0]][d];
}

} // luma
