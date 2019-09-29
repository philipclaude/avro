// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "common/error.h"

#include "mesh/vertices.h"

//#include "numerics/densmat.h"
//#include "numerics/determinant.h"
//#include "numerics/expansion.h"
//#include "numerics/functions.h"
#include "numerics/geometry.h"
//#include "numerics/orient4d_filter.h"

//#include <triangle/predicates.h>

//typedef ursa::real_t REAL;
//#include <tetgen1.5.0/predicates.h>

namespace ursa
{

namespace numerics
{

/*
real_t
orient4d( const real_t* p0 , const real_t* p1 , const real_t* p2 , const real_t* p3 , const real_t* p4 )
{
  using namespace GEO;
  using namespace GEO::PCK;

  int s = orient4d_filter(p0,p1,p2,p3,p4);
  if (s!=FPG_UNCERTAIN_VALUE)
  {
    return orient4dfast(p0,p1,p2,p3,p4);
  }
  const expansion& a11 = expansion_diff( p1[0] , p0[0] );
  const expansion& a12 = expansion_diff( p2[0] , p0[0] );
  const expansion& a13 = expansion_diff( p3[0] , p0[0] );
  const expansion& a14 = expansion_diff( p4[0] , p0[0] );

  const expansion& a21 = expansion_diff( p1[1] , p0[1] );
  const expansion& a22 = expansion_diff( p2[1] , p0[1] );
  const expansion& a23 = expansion_diff( p3[1] , p0[1] );
  const expansion& a24 = expansion_diff( p4[1] , p0[1] );

  const expansion& a31 = expansion_diff( p1[2] , p0[2] );
  const expansion& a32 = expansion_diff( p2[2] , p0[2] );
  const expansion& a33 = expansion_diff( p3[2] , p0[2] );
  const expansion& a34 = expansion_diff( p4[2] , p0[2] );

  const expansion& a41 = expansion_diff( p1[3] , p0[3] );
  const expansion& a42 = expansion_diff( p2[3] , p0[3] );
  const expansion& a43 = expansion_diff( p3[3] , p0[3] );
  const expansion& a44 = expansion_diff( p4[3] , p0[3] );

  const expansion& Delta = expansion_det4x4( a11 , a12 , a13 , a14 , a21 , a22 , a23 , a24 , a31 ,a32,a33,a34 , a41,a42,a43,a44 );
  if (Delta.sign()==GEO::ZERO) return 0.0;
  return Delta.value();
}
*/

real_t
distance2( const real_t* x , const real_t* y , const coord_t dim )
{
  real_t d = 0.;
  for (coord_t i=0;i<dim;i++)
    d += ( x[i] -y[i] )*( x[i] -y[i] );
  return d;
}

real_t
distance( const real_t* x , const real_t* y , const coord_t dim )
{
  return sqrt(distance2(x,y,dim));
}
/*

real_t
volume2( const std::vector<real_t*>& x )
{
  return .5*orient2d(x[0],x[1],x[2]);
}

real_t
volume3( const std::vector<real_t*>& x )
{
  return -orient3d(x[0],x[1],x[2],x[3])/6.;
}

real_t
volume4( const std::vector<real_t*>& x )
{
  real_t dv = orient4d(x[0],x[1],x[2],x[3],x[4])/24.;
  return dv;
}

real_t
volume( const std::vector<real_t*>& x , const coord_t dim )
{
  coord_t n = x.size() -1;
  if (n==2 && dim==2) return fabs(volume2(x));
  if (n==3 && dim==3) return fabs(volume3(x));
  if (n==4 && dim==4) return fabs(volume4(x));
  std::vector<const real_t*> y(x.begin(),x.end()); // use the function below
  return volume( y , dim );
  avro_assert_not_reached;
  return 0.;
}

real_t
signedVolume( const std::vector<real_t*>& x , const coord_t dim )
{
  coord_t n = x.size() -1;
  if (n==2 && dim==2) return volume2(x);
  if (n==3 && dim==3) return volume3(x);
  if (n==4 && dim==4) return volume4(x);
  avro_assert_not_reached;
  return 0.;
}

real_t
volume( const std::vector<const real_t*>& x , const coord_t dim )
{
  // assume this is a simplex because there's no other way to compute it
  coord_t n = x.size() -1;

  numerics::densMat<real_t> B(n+2,n+2);

  B(0,0) = 0;
  for (index_t i=1;i<index_t(n+2);i++)
  {
    B(i,0) = 1.;
    B(0,i) = 1.;
  }

  for (index_t i=0;i<index_t(n+1);i++)
  for (index_t j=0;j<index_t(n+1);j++)
    B(i+1,j+1) = distance2(x[i],x[j],dim);

  real_t f = pow(-1.,n+1)/pow(2.,n);
  f /= pow( real_t(numerics::factorial(n)),2. );
  real_t v = sqrt(fabs(f*numerics::determinant(B)));
  return v;
}

void
barycentric( const real_t* p , const std::vector<const real_t*>& x , const coord_t dim , std::vector<real_t>& alpha )
{
  index_t n = x.size() -1;
  alpha.resize(n+1);

  real_t V = volume(x,dim);

  for (index_t k=0;k<n+1;k++)
  {
    std::vector<const real_t*> f = x;
    f.erase( f.begin() + k );

    f.push_back( p );

    real_t vk = volume(f,dim);

    alpha[k] = vk/V;
  }
}

void
barycentric_signed( real_t* p , const std::vector<real_t*>& x , const coord_t dim , std::vector<real_t>& alpha )
{
  index_t n = x.size() -1;
  alpha.resize(n+1);

  real_t V = signedVolume(x,dim);
  if (V<=0.0) printf("warning: V = %1.12e\n",V);
  //avro_assert( V>0.0 );

  for (index_t k=0;k<n+1;k++)
  {
    std::vector<real_t*> f = x;
    f[k] = p;

    real_t vk = signedVolume(f,dim);

    alpha[k] = vk/V;
  }

}

void
centroid( const index_t* v0 , const index_t nv , const Vertices& vertices , std::vector<real_t>& xc )
{
  avro_assert_msg( vertices.dim()==xc.size() ,
    "dim(v) = %u, dim(xc) = %lu" , vertices.dim() , xc.size() );

  std::fill( xc.begin() , xc.end() , 0. );
  for (index_t j=0;j<nv;j++)
  for (coord_t d=0;d<vertices.dim();d++)
    xc[d] += vertices(*(v0+j),d)/nv;
}

void
normal( const std::vector<real_t*>& x , real_t* n , const coord_t dim )
{
  // we want a single vector orthogonal to the simplex in x
  // in a higher dimensional space, we would have multiple vectors
  avro_assert_msg( x.size()==dim , "|x| = %lu, dim = %u" , x.size() , dim );

  numerics::densMat<double> A(dim-1,dim);
  A.zeros();

  for (index_t k=0;k<index_t(dim-1);k++)
  for (index_t j=0;j<dim;j++)
    A(k,j) = x[k+1][j] -x[0][j];

  numerics::densMat<double> K(dim-1,1);
  int info = A.kernel(K);
  if (info!=0) A.display();
  avro_assert_msg( info>=0 , " info = %d " , info );

  for (index_t j=0;j<dim;j++)
    n[j] = K(j,0);

  normalize(n,dim);
}


void
orthogonal( const std::vector<real_t*>& x , real_t* n , const coord_t dim )
{
  // we want a single vector orthogonal to the vector in x
  // in a higher dimensional space, we would have multiple vectors
  avro_assert_msg( x.size()==index_t(dim-1) ,
                    "|x| = %lu, dim = %u" , x.size() , dim );

  numerics::densMat<double> A(dim-1,dim);
  A.zeros();

  for (index_t k=0;k<index_t(dim-1);k++)
  for (index_t j=0;j<dim;j++)
    A(k,j) = x[k][j];

  numerics::densMat<double> K(dim-1,1);
  int info = A.kernel(K);
  if (info!=0) A.display();
  avro_assert_msg( info>=0 , " info = %d " , info );

  for (index_t j=0;j<dim;j++)
    n[j] = K(j,0);

  normalize(n,dim);
}

void
orientNormal( const std::vector<real_t*>& x , real_t *n , const coord_t dim )
{
  // create a new point using the first ppoint and the normal
  std::vector<real_t> xp( dim , 0. );
  for (coord_t d=0;d<dim;d++)
    xp[d] = x[0][d] +n[d];

  // check if the orientation is correct
  bool flip = false;
  if (dim==2)
  {
    if (orient2d(x[0],x[1],xp.data())<0)
      flip = true;
  }
  else if (dim==3)
  {
    if (orient3d(x[0],x[1],x[2],xp.data())<0)
      flip = true;
  }
  else
  {
    printf("dimension = %d not implemented\n",int(dim));
    avro_assert_not_reached;
  }

  // flip the normal if necessary
  if (flip)
  {
    for (coord_t d=0;d<dim;d++)
      n[d] = -n[d];
  }

}

void
normalize( real_t* x , const coord_t dim )
{
  real_t l = 0.;
  for (coord_t d=0;d<dim;d++)
    l += x[d]*x[d];
  l = std::sqrt(l);
  for (coord_t d=0;d<dim;d++)
    x[d] /= l;
}

void
axpb( const real_t a , const real_t* x , const real_t* b , const coord_t dim , real_t* y )
{
  for (coord_t d=0;d<dim;d++)
    y[d] = a*x[d] +b[d];
}

void
vector( const real_t* x0 , const real_t* x1 , const coord_t dim , real_t* v )
{
  for (coord_t d=0;d<dim;d++)
    v[d] = x1[d] -x0[d];
}
*/

} // numerics

} // ursa
