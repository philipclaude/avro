// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(ELEMENTARYREFLECTOR_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include <limits>
#include <cmath> // fabs

#include "ElementaryReflector.h"

#include "MatrixD.h"
#include "VectorD.h"
#include "MatrixD_Transpose.h"

#include "tinymat/dense/tools/dot.h"
#include "tinymat/dense/tools/norm2.h"

namespace tinymat 
{
namespace DLA
{

#define SIGN( x, y ) (y) > 0 ? (x) : -(x)

//=====================================================================
template< class T, class T2 >
void ElementaryReflector_impl<T,T2>::apply( const VectorDView< T >& v, const T& tau, MatrixDView< T2 >& C )
/*
Based on the LAPACK routine DLARF

ApplyElementaryReflector applies a real elementary reflector H to a real m by n matrix
C, from either the left or the right. H is represented in the form

      H = I - tau * v * v**T

where tau is a real scalar and v is a real vector.

If tau = 0, then H is taken to be the unit matrix.
*/
{

  int lastv = -1;
  int lastc = -1;
  if ( tau != 0 )
  {
    // Set up variables for scanning V.  lastv begins pointing to the end of V.
    lastv = v.m();

    // Look for the last non-zero row in V.
    while ( lastv > 0 && v[lastv-1] == 0 )
      lastv--;

    // Scan for the last non-zero column in C(1:lastv,:).
    lastc = C.n();
    //well... this is probably more general than we care for...
  }

  // Form  H * C

  if ( lastv >= 0 )
  {
    MatrixDView< T2 > c = C.sub(0,0,lastv,lastc);
    MatrixD<T2> w(1,lastc);
    MatrixDView< T > vsub = v.sub(0,0,lastv,1);

    //
    // w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1)
    //
    w = Transpose(vsub)*c;
    //
    // C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**T
    //
    c -= tau*vsub*w;
  }
}

//Compute sqrt( x'*x ) following the algorithm in LAPACK DNRM2
template< class T >
T dnrm2(const VectorDView< T >& x)
{
  T scale = 0;
  T ssq = 0;
  const int m = x.m();
  for ( int i = 0; i < m; i++ )
    if ( x[i] != 0 )
    {
      T absxi = fabs(x[i]);
      if ( scale < absxi )
      {
        ssq = T(1) + ssq* (scale/absxi)*(scale/absxi);
        scale = absxi;
      }
      else
        ssq += (absxi/scale)*(absxi/scale);
    }

  return scale*sqrt(ssq);
}


template< class T, class T2 >
void ElementaryReflector_impl<T,T2>::compute( T& alpha, VectorDView< T >& x, T& tau )
/*
Based on the LAPACK routine DLARFG

ElementaryReflector generates a real elementary reflector H of order n, such
that

      H * ( alpha ) = ( beta ),   H**T * H = I.
          (   x   )   (   0  )

where alpha and beta are scalars, and x is an (n-1)-element real
vector. H is represented in the form

      H = I - tau * ( 1 ) * ( 1 v**T ) ,
                    ( v )

where tau is a real scalar and v is a real (n-1)-element
vector.

If the elements of x are all zero, then tau = 0 and H is taken to be
the unit matrix.

Otherwise  1 <= tau <= 2.
*/
{

  //Compute sqrt( x'*x )
  T xnorm = dnrm2(x);

  if ( xnorm == 0 )
  {
    //H  =  I
    tau = 0;
  }
  else
  {
    //general case
    T beta = -SIGN( norm2(alpha,xnorm), alpha );
    T safmin = std::numeric_limits<T>::min();
    T small = T(1) / std::numeric_limits<T>::max();
    if ( small >= safmin )
    {
      //Use SMALL plus a bit, to avoid the possibility of rounding
      //causing overflow when computing  1/sfmin.
      safmin = small*( T(1)+std::numeric_limits<T>::epsilon() );
    }

    int knt = 0;
    if ( fabs( beta ) < safmin )
    {
      //xnorm, beta may be inaccurate; scale X and recompute them

      while ( fabs( beta ) < safmin )
      {
        T rsafmn = T(1) / safmin;
        knt = knt + 1;
        x *= rsafmn;
        beta = beta*rsafmn;
        alpha = alpha*rsafmn;
      }
      //New beta is at most 1, at least safmin
      xnorm = dnrm2(x);
      beta = -SIGN( norm2(alpha,xnorm), alpha );
    }

    tau = ( beta-alpha ) / beta;
    x /= ( alpha-beta );

    //If alpha is subnormal, it may lose relative accuracy
    for ( int j = 0; j < knt; j++)
      beta *= safmin;

    alpha = beta;
  }
}


template< int M, int N, class T1, class T2 >
struct ElementaryReflector_impl< MatrixS<M,N,T1>, T2 >
{
  typedef MatrixS<M,N,T1> T;
  static void compute( T& alpha, VectorDView< T >& x, T& tau )
  {
    SANS_ASSERT(false); //Yeah I have not figured out how to do this for matrices yet
    //Compute sqrt( x'*x )
    //T xnorm = Transpose(x)*x;
    //T beta = -SIGN( sqrt(alpha*alpha + xnorm), alpha );

    //tau = InverseQR::Inverse( beta )*( beta-alpha );
    //x = InverseQR::Inverse( alpha-beta )*x;
  }

  static void apply( const VectorDView< T >& v, const T& tau, MatrixDView< T2 >& C )
  {
    SANS_ASSERT(false); //Yeah I have not figured out how to do this for matrices yet

  }
};

} //namespace DLA
} //namespace tinymat 
