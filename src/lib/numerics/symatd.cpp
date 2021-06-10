// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "numerics/mat.h"
#include "numerics/determinant.h"
#include "numerics/dual.h"
#include "numerics/symatd.h"

#include "numerics/surreal/config.h"
#include "numerics/surreal/SurrealD.h"
#include "numerics/surreal/SurrealS.h"

#include <stdio.h>
#include <math.h>

namespace avro
{

namespace numerics
{

template<typename type>
symd<type>::symd( const coord_t _n ) :
	n_(_n)
{
	// constructor from size
	data_.resize( nb() );
	std::fill( data_.begin() , data_.end() , 0. );
}

template<typename type>
symd<type>::symd( const std::vector<type>& _data ) :
  data_(_data)
{
  // constructor from matrix data itself
  index_t m = data_.size();
  n_ = ( -1 +::sqrt( 8*m +1 ) )/2;
}

template<typename type>
symd<type>::symd( const vecd<type>& lambda , const matd<type>& q ) :
  n_( lambda.m() )
{
  // constructor from eigendecomposition
  data_.resize( nb() );
  //std::fill( data_.begin() , data_.end() , 0. );

	(*this) = q * diag(lambda) * transpose(q);

	/*
	// compute L*q'
  matd<type> lqt;
  lqt.copy(q,true); // true for transpose
  for (index_t k=0;k<n_;k++)
  for (index_t j=0;j<n_;j++)
    lqt(k,j) *= lambda(k);

  // compute q * L * q'
  matd<type> S(n_,n_);
  q.multiply(lqt,S);

  // save the result
  for (index_t i=0;i<n_;i++)
  for (index_t j=0;j<n_;j++)
    operator()(i,j) = S(i,j);
	*/
}

template<typename type>
symd<type>::symd( const std::pair< vecd<type> , matd<type> >& decomp ) :
  symd( decomp.first , decomp.second )
{
}

template<typename type>
symd<type>::symd( const matd<type>& M ) :
  n_(M.n())
{
  avro_assert( M.m() == M.n() );
  data_.resize(nb());
  for (index_t i=0;i<n_;i++)
  for (index_t j=0;j<n_;j++)
    operator()(i,j) = M(i,j);
}

template<typename type>
template<typename typeR>
symd<type>::symd( const symd<typeR>& A ) :
	n_(A.n())
{
	data_.resize(nb());
	for (index_t i=0;i<n_;i++)
	for (index_t j=0;j<n_;j++)
		operator()(i,j) = type(A(i,j));
}

#if USE_SURREAL
template symd<SurrealS<1>>::symd( const symd<real_t>& A );
template symd<SurrealS<3>>::symd( const symd<real_t>& A );
template symd<SurrealS<6>>::symd( const symd<real_t>& A );
template symd<SurrealS<10>>::symd( const symd<real_t>& A );
#endif

template<typename type>
void
symd<type>::fromeig(const vecd<type>& lambda , const matd<type>& q )
{
  n_ = lambda.m();

  // constructor from eigendecomposition
  data_.resize( nb() );
  std::fill( data_.begin() , data_.end() , 0. );

	(*this) = q * diag(lambda) * transpose(q);

	/*
  // compute L*q'
  matd<type> lqt;
  lqt.copy(q,true); // true for transpose
  for (index_t k=0;k<n_;k++)
  for (index_t j=0;j<n_;j++)
    lqt(k,j) *= lambda[k];

  // compute q * L * q'
  matd<type> S(n_,n_);
  q.multiply(lqt,S);

  // save the result
  for (index_t i=0;i<n_;i++)
  for (index_t j=0;j<n_;j++)
    operator()(i,j) = S(i,j);
	*/
}

template<typename type>
type
symd<type>::quadratic_form( const type* x ) const
{
  // computes the quadratic form x'*(this)*x

  // y = (this)*x
  std::vector<type> y(n_);
  for (coord_t i=0;i<n_;i++) // loop through rows of (this)
  {
    y[i] = 0.;
    for (coord_t j=0;j<n_;j++) // loop through colums of (this)
    {
      y[i] += operator()(i,j)*x[j];
    }
  }

  // x' * y
  type f = 0.;
  for (coord_t i=0;i<n_;i++)
  {
    f += x[i]*y[i];
  }
  return f;
}

/*
template<typename type>
type
symd<type>::quadraticFormReal( real_t* x ) const
{
  // computes the quadratic form x'*(this)*x

  // y = (this)*x
  std::vector<type> y(n_);
  for (coord_t i=0;i<n_;i++) // loop through rows of (this)
  {
    y[i] = 0.;
    for (coord_t j=0;j<n_;j++) // loop through colums of (this)
    {
      y[i] += operator()(i,j)*x[j];
    }
  }

  // x' * y
  type f = 0.;
  for (coord_t i=0;i<n_;i++)
  {
    f += x[i]*y[i];
  }
  return f;
}
*/

template<typename type>
type
symd<type>::determinant() const
{
  if (n_==2)
		return data_[0]*data_[2] -data_[1]*data_[1];
  if (n_==3)
  {
		return  operator()(0,0)*(operator()(2,2)*operator()(1,1)-operator()(2,1)*operator()(1,2))
	         -operator()(1,0)*(operator()(2,2)*operator()(0,1)-operator()(2,1)*operator()(0,2))
	         +operator()(2,0)*(operator()(1,2)*operator()(0,1)-operator()(1,1)*operator()(0,2));
  }
	assert( n_ == 4 );
	type X1_1=operator()(0,0),X1_2=operator()(0,1),X1_3=operator()(0,2),X1_4=operator()(0,3);
	type X2_1=operator()(1,0),X2_2=operator()(1,1),X2_3=operator()(1,2),X2_4=operator()(1,3);
	type X3_1=operator()(2,0),X3_2=operator()(2,1),X3_3=operator()(2,2),X3_4=operator()(2,3);
	type X4_1=operator()(3,0),X4_2=operator()(3,1),X4_3=operator()(3,2),X4_4=operator()(3,3);
	return
        X1_1*X2_2*X3_3*X4_4 - X1_1*X2_2*X3_4*X4_3 - X1_1*X2_3*X3_2*X4_4
      + X1_1*X2_3*X3_4*X4_2 + X1_1*X2_4*X3_2*X4_3 - X1_1*X2_4*X3_3*X4_2
      - X1_2*X2_1*X3_3*X4_4 + X1_2*X2_1*X3_4*X4_3 + X1_2*X2_3*X3_1*X4_4
      - X1_2*X2_3*X3_4*X4_1 - X1_2*X2_4*X3_1*X4_3 + X1_2*X2_4*X3_3*X4_1
      + X1_3*X2_1*X3_2*X4_4 - X1_3*X2_1*X3_4*X4_2 - X1_3*X2_2*X3_1*X4_4
      + X1_3*X2_2*X3_4*X4_1 + X1_3*X2_4*X3_1*X4_2 - X1_3*X2_4*X3_2*X4_1
      - X1_4*X2_1*X3_2*X4_3 + X1_4*X2_1*X3_3*X4_2 + X1_4*X2_2*X3_1*X4_3
      - X1_4*X2_2*X3_3*X4_1 - X1_4*X2_3*X3_1*X4_2 + X1_4*X2_3*X3_2*X4_1;
}

/*
template<typename type>
symd<type>
symd<type>::sandwich( const symd<type>& B ) const
{
	symd<type> C(n_);

	avro_assert( n_>=2 && n_<=4 );

	type A1_1 = operator()(0,0); type A1_2 = operator()(0,1);
	type A2_1 = operator()(1,0); type A2_2 = operator()(1,1);
	type B1_1 = B(0,0); type B1_2 = B(0,1);
	type B2_1 = B(1,0); type B2_2 = B(1,1);

	if (n_==2)
	{
		C(0,0) = B1_1*(A1_1*B1_1+A2_1*B1_2)+B2_1*(A1_2*B1_1+A2_2*B1_2);
		C(0,1) = B1_2*(A1_1*B1_1+A2_1*B1_2)+B2_2*(A1_2*B1_1+A2_2*B1_2);
		C(1,1) = B1_2*(A1_1*B2_1+A2_1*B2_2)+B2_2*(A1_2*B2_1+A2_2*B2_2);
	}
	else if (n_==3)
	{
		type A1_3 = operator()(0,2); type A2_3 = operator()(1,2); type A3_3 = operator()(2,2);
		type A3_1 = A1_3; type A3_2 = A2_3;

		type B1_3 = B(0,2); type B2_3 = B(1,2); type B3_3 = B(2,2);
		type B3_1 = B1_3; type B3_2 = B2_3;

		C(0,0) = B1_1*(A1_1*B1_1+A2_1*B1_2+A3_1*B1_3)+B2_1*(A1_2*B1_1+A2_2*B1_2+A3_2*B1_3)+B3_1*(A1_3*B1_1+A2_3*B1_2+A3_3*B1_3);
    C(0,1) = B1_2*(A1_1*B1_1+A2_1*B1_2+A3_1*B1_3)+B2_2*(A1_2*B1_1+A2_2*B1_2+A3_2*B1_3)+B3_2*(A1_3*B1_1+A2_3*B1_2+A3_3*B1_3);
    C(0,2) = B1_3*(A1_1*B1_1+A2_1*B1_2+A3_1*B1_3)+B2_3*(A1_2*B1_1+A2_2*B1_2+A3_2*B1_3)+B3_3*(A1_3*B1_1+A2_3*B1_2+A3_3*B1_3);
    C(1,1) = B1_2*(A1_1*B2_1+A2_1*B2_2+A3_1*B2_3)+B2_2*(A1_2*B2_1+A2_2*B2_2+A3_2*B2_3)+B3_2*(A1_3*B2_1+A2_3*B2_2+A3_3*B2_3);
    C(1,2) = B1_3*(A1_1*B2_1+A2_1*B2_2+A3_1*B2_3)+B2_3*(A1_2*B2_1+A2_2*B2_2+A3_2*B2_3)+B3_3*(A1_3*B2_1+A2_3*B2_2+A3_3*B2_3);
    C(2,2) = B1_3*(A1_1*B3_1+A2_1*B3_2+A3_1*B3_3)+B2_3*(A1_2*B3_1+A2_2*B3_2+A3_2*B3_3)+B3_3*(A1_3*B3_1+A2_3*B3_2+A3_3*B3_3);
	}
	else if (n_==4)
	{
		type A1_3 = operator()(0,2); type A2_3 = operator()(1,2); type A3_3 = operator()(2,2);
		type A3_1 = A1_3; type A3_2 = A2_3;

		type B1_3 = B(0,2); type B2_3 = B(1,2); type B3_3 = B(2,2);
		type B3_1 = B1_3; type B3_2 = B2_3;

		type A1_4 = operator()(0,3); type A2_4 = operator()(1,3); type A3_4 = operator()(2,3); type A4_4 = operator()(3,3);
		type A4_1 = A1_4; type A4_2 = A2_4; type A4_3 = A3_4;

		type B1_4 = B(0,3); type B2_4 = B(1,3); type B3_4 = B(2,3); type B4_4 = B(3,3);
		type B4_1 = B1_4; type B4_2 = B2_4; type B4_3 = B3_4;

		C(0,0) = B1_1*(A1_1*B1_1+A2_1*B1_2+A3_1*B1_3+A4_1*B1_4)+B2_1*(A1_2*B1_1+A2_2*B1_2+A3_2*B1_3+A4_2*B1_4)+B3_1*(A1_3*B1_1+A2_3*B1_2+A3_3*B1_3+A4_3*B1_4)+B4_1*(A1_4*B1_1+A2_4*B1_2+A3_4*B1_3+A4_4*B1_4);
    C(0,1) = B1_2*(A1_1*B1_1+A2_1*B1_2+A3_1*B1_3+A4_1*B1_4)+B2_2*(A1_2*B1_1+A2_2*B1_2+A3_2*B1_3+A4_2*B1_4)+B3_2*(A1_3*B1_1+A2_3*B1_2+A3_3*B1_3+A4_3*B1_4)+B4_2*(A1_4*B1_1+A2_4*B1_2+A3_4*B1_3+A4_4*B1_4);
    C(0,2) = B1_3*(A1_1*B1_1+A2_1*B1_2+A3_1*B1_3+A4_1*B1_4)+B2_3*(A1_2*B1_1+A2_2*B1_2+A3_2*B1_3+A4_2*B1_4)+B3_3*(A1_3*B1_1+A2_3*B1_2+A3_3*B1_3+A4_3*B1_4)+B4_3*(A1_4*B1_1+A2_4*B1_2+A3_4*B1_3+A4_4*B1_4);
    C(0,3) = B1_4*(A1_1*B1_1+A2_1*B1_2+A3_1*B1_3+A4_1*B1_4)+B2_4*(A1_2*B1_1+A2_2*B1_2+A3_2*B1_3+A4_2*B1_4)+B3_4*(A1_3*B1_1+A2_3*B1_2+A3_3*B1_3+A4_3*B1_4)+B4_4*(A1_4*B1_1+A2_4*B1_2+A3_4*B1_3+A4_4*B1_4);
    C(1,1) = B1_2*(A1_1*B2_1+A2_1*B2_2+A3_1*B2_3+A4_1*B2_4)+B2_2*(A1_2*B2_1+A2_2*B2_2+A3_2*B2_3+A4_2*B2_4)+B3_2*(A1_3*B2_1+A2_3*B2_2+A3_3*B2_3+A4_3*B2_4)+B4_2*(A1_4*B2_1+A2_4*B2_2+A3_4*B2_3+A4_4*B2_4);
    C(1,2) = B1_3*(A1_1*B2_1+A2_1*B2_2+A3_1*B2_3+A4_1*B2_4)+B2_3*(A1_2*B2_1+A2_2*B2_2+A3_2*B2_3+A4_2*B2_4)+B3_3*(A1_3*B2_1+A2_3*B2_2+A3_3*B2_3+A4_3*B2_4)+B4_3*(A1_4*B2_1+A2_4*B2_2+A3_4*B2_3+A4_4*B2_4);
    C(1,3) = B1_4*(A1_1*B2_1+A2_1*B2_2+A3_1*B2_3+A4_1*B2_4)+B2_4*(A1_2*B2_1+A2_2*B2_2+A3_2*B2_3+A4_2*B2_4)+B3_4*(A1_3*B2_1+A2_3*B2_2+A3_3*B2_3+A4_3*B2_4)+B4_4*(A1_4*B2_1+A2_4*B2_2+A3_4*B2_3+A4_4*B2_4);
    C(2,2) = B1_3*(A1_1*B3_1+A2_1*B3_2+A3_1*B3_3+A4_1*B3_4)+B2_3*(A1_2*B3_1+A2_2*B3_2+A3_2*B3_3+A4_2*B3_4)+B3_3*(A1_3*B3_1+A2_3*B3_2+A3_3*B3_3+A4_3*B3_4)+B4_3*(A1_4*B3_1+A2_4*B3_2+A3_4*B3_3+A4_4*B3_4);
    C(2,3) = B1_4*(A1_1*B3_1+A2_1*B3_2+A3_1*B3_3+A4_1*B3_4)+B2_4*(A1_2*B3_1+A2_2*B3_2+A3_2*B3_3+A4_2*B3_4)+B3_4*(A1_3*B3_1+A2_3*B3_2+A3_3*B3_3+A4_3*B3_4)+B4_4*(A1_4*B3_1+A2_4*B3_2+A3_4*B3_3+A4_4*B3_4);
    C(3,3) = B1_4*(A1_1*B4_1+A2_1*B4_2+A3_1*B4_3+A4_1*B4_4)+B2_4*(A1_2*B4_1+A2_2*B4_2+A3_2*B4_3+A4_2*B4_4)+B3_4*(A1_3*B4_1+A2_3*B4_2+A3_3*B4_3+A4_3*B4_4)+B4_4*(A1_4*B4_1+A2_4*B4_2+A3_4*B4_3+A4_4*B4_4);
	}
	else
		avro_implement;
	return C;
}
*/

template<typename type>
symd<type>
symd<type>::exp() const
{
	std::pair< vecd<type>,matd<type> > decomp = eig();
	for (index_t k=0;k<n_;k++)
		decomp.first(k) = ::exp( decomp.first(k) );
	return symd<type>(decomp);
}

template<typename type>
symd<type>
symd<type>::log() const
{
	std::pair< vecd<type>,matd<type> > decomp = eig();
	for (index_t k=0;k<n_;k++)
		decomp.first(k) = ::log( decomp.first(k) );
	return symd<type>(decomp);
}

template<typename type>
symd<type>
symd<type>::pow( real_t p ) const
{
	std::pair< vecd<type>,matd<type> > decomp = eig();
	for (index_t k=0;k<n_;k++)
		decomp.first(k) = ::pow( decomp.first(k) , p );
	return symd<type>(decomp);
}

template<typename type>
symd<type>
symd<type>::sqrt() const
{
	std::pair< vecd<type>,matd<type> > decomp = eig();
	for (index_t k=0;k<n_;k++)
		decomp.first(k) = ::sqrt( decomp.first(k) );
	return symd<type>(decomp);
}

// eigenvalues and eigenvectors
template<typename type>
std::pair< vecd<type>,matd<type> >
symd<type>::eig() const
{
  if (n_==2) return __eigivens__();
  if (n_==3 || n_==4) return __eigivens__();
  printf("unknown eigendecomposition for %ux%u tensors\n",n_,n_);
  avro_assert_not_reached;
	return { vecd<type>(n_) , matd<type>(n_,n_) };
}

template<typename type>
std::pair< vecd<type>,matd<type> >
symd<type>::__eig2__() const
{
  vecd<type> d(n_);
  matd<type> Q(n_,n_);
  type sqDelta,dd,trm,vnorm;
  dd  = data_[0] -data_[2];
  trm = data_[0] +data_[2];
  sqDelta = ::sqrt( dd*dd +4.0*data_[1]*data_[1]);

  d[0] = 0.5*(trm -sqDelta);

  real_t tol = 1e-12;
  if (sqDelta < tol )
  {
    d[1] = d[0];
    Q(0,0) = 1.;
    Q(0,1) = 0.;
    Q(1,0) = 0.;
    Q(1,1) = 1.;
    return std::make_pair(d,Q);
  }

  Q(0,0) = data_[1];
  Q(0,1) = (d[0] -data_[0]);
  vnorm = ::sqrt( Q(0,0)*Q(0,0) +Q(0,1)*Q(0,1) );

  if (vnorm<tol)
  {
    Q(0,0) = (d[0] -data_[2]);
    Q(0,1) = data_[1];
    vnorm = ::sqrt( Q(0,0)*Q(0,0) +Q(0,1)*Q(0,1) );
  }

  avro_assert( vnorm>tol );

  vnorm = 1./vnorm;
  Q(0,0) *= vnorm;
  Q(0,1) *= vnorm;

  Q(1,0) = -Q(0,1);
  Q(1,1) = Q(0,0);

  d[1] =    data_[0]*Q(1,0)*Q(1,0)
        +2.*data_[1]*Q(1,0)*Q(1,1)
        +   data_[2]*Q(1,1)*Q(1,1);

  return std::make_pair(d,Q);
}

template<typename type>
std::pair< vecd<type>,matd<type> >
symd<type>::__eigivens__() const
{
  vecd<type> L(n_);
  matd<type> E(n_,n_);

  type sd,so;
  type s,c,t;
  type g,h,z,theta;
  type thresh;

  E.eye();

  symd A(n_);
  A.copy(*this);

  for (coord_t i=0;i<n_;i++)
    L[i] = A(i,i);

  // calculate sq(tr(A))
  sd = 0.;
  for (coord_t i=0;i<n_;i++)
    sd += fabs(L[i]);
  sd = sd*sd;

  for (index_t iter=0;iter<50;iter++)
  {
    // test for convergence
    so = 0.;
    for (coord_t p=0;p<n_;p++)
      for (coord_t q=p+1;q<n_;q++)
        so += fabs(A(p,q));

    if (so==0.)
      return std::make_pair(L,E);

    if (iter<4)
      thresh = 0.2*so/type(n_*n_);
    else
      thresh = 0.;

    // sweep
    for (coord_t p=0;p<n_;p++)
    {
      for (int q=p+1;q<n_;q++)
      {
        g = 100.*fabs(A(p,q));
        if (iter>4 && fabs(L[p])+g == fabs(L[p]) && fabs(L[q])+g==fabs(L[q]))
          A(p,q) = 0.;
        else if (fabs(A(p,q)) > thresh)
        {
          // calculate Jacobi transformation
          h = L[q] -L[p];
          if (fabs(h)+g==fabs(h))
          {
            t = A(p,q)/h;
          }
          else
          {
            theta = 0.5*h/A(p,q);
            if (theta<0.)
              t = -1./( ::sqrt(1. +theta*theta) -theta );
            else
              t = 1./( ::sqrt(1. +theta*theta) +theta );
          }
          c = 1./::sqrt(1. +t*t);
          s = t*c;
          z = t*A(p,q);

          // apply Jacobi transformation
          A(p,q) = 0.;
          L[p] -= z;
          L[q] += z;
          for (coord_t r=0;r<p;r++)
          {
            t = A(r,p);
            A(r,p) = c*t -s*A(r,q);
            A(r,q) = s*t +c*A(r,q);
          }

          for (coord_t r=p+1;r<q;r++)
          {
            t = A(p,r);
            A(p,r) = c*t -s*A(r,q);
            A(r,q) = s*t +c*A(r,q);
          }

          for (coord_t r=q+1;r<n_;r++)
          {
            t = A(p,r);
            A(p,r) = c*t -s*A(q,r);
            A(q,r) = s*t +c*A(q,r);
          }

          // update eigenvectors
          for (coord_t r=0;r<n_;r++)
          {
            t = E(r,p);
            E(r,p) = c*t -s*E(r,q);
            E(r,q) = s*t +c*E(r,q);
          }
        }
      }
    }
  } // iteration loop

	printf("givens rotation failed :(\n");
	display();

	avro_assert_not_reached;

  return std::make_pair(L,E);
}

template<typename type>
bool
symd<type>::check()
{
  // check if the tensor is indeed symd be ensuring the eigenvalues are > 0
  std::pair< vecd<type>,matd<type> > E = eig();
  for (index_t k=0;k<E.first.m();k++)
  {
    if (E.first[k]<1e-16)
      return false;
  }
  return true;
}

template<typename type>
void
symd<type>::display( const std::string& title ) const
{
	avro_assert_not_reached;
}

template<typename type>
void
symd<type>::matlabize( const std::string& title ) const
{
	avro_assert_not_reached;
}

template<typename type>
void
symd<type>::forunit( const std::string& title ) const
{
	avro_assert_not_reached;
}

template<>
void
symd<real_t>::display( const std::string& title ) const
{
  if (!title.empty()) printf("%s\n",title.c_str());
  else printf("SPT:\n");
  for (index_t i=0;i<n_;i++)
  {
    printf("[ ");
    for (index_t j=0;j<n_;j++)
      printf("%.16e ",operator()(i,j));
    printf("]\n");
  }
}

template<>
void
symd<real_t>::matlabize( const std::string& title ) const
{
	if (!title.empty()) printf("%s = [",title.c_str());
  else printf("A = [");
  for (index_t i=0;i<n_;i++)
  {
    for (index_t j=0;j<n_;j++)
		{
      printf("%.16e",operator()(i,j));
			if (int(j)<n_-1) printf(",");
			else
			{
				if (int(i)<n_-1)
					printf(";");
			}
		}
  }
	printf("];\n");
}

template<>
void
symd<real_t>::forunit( const std::string& title ) const
{
	avro_assert_msg( !title.empty() , "provide a variable name!" );
	for (index_t i=0;i<n_;i++)
	for (index_t j=0;j<=i;j++)
		printf("%s(%lu,%lu) = %.16e;\n",title.c_str(),i,j,operator()(i,j));
}

// framework-specific functions
template<typename type>
void
LogEuclidean<type>::interpolate( const std::vector<type>& alpha ,
                         const std::vector<symd<type>>& tensors , symd<type>& T )
{
  avro_assert( alpha.size()==tensors.size() );

  T.zero();
  for (index_t k=0;k<tensors.size();k++)
	{
    T = T +tensors[k].log()*alpha[k];
	}
  T = T.exp();
}

template class symd<real_t>;
template class LogEuclidean<real_t>;

template class symd<dual>;
template class LogEuclidean<dual>;

#if 0 //USE_SURREAL
template class symd<SurrealS<1>>;
template class LogEuclidean<SurrealS<1>>;

template class symd<SurrealS<3>>;
template class LogEuclidean<SurrealS<3>>;

template class symd<SurrealS<6>>;
template class LogEuclidean<SurrealS<6>>;

template class symd<SurrealS<10>>;
template class LogEuclidean<SurrealS<10>>;

#endif

} // numerics

} // avro
