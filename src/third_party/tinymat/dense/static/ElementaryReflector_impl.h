// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(MATRIXS_ELEMENTARYREFLECTOR_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include <limits>
#include <cmath> // fabs

#include "ElementaryReflector.h"

#include "MatrixS.h"
#include "VectorS.h"
#include "MatrixS_Transpose.h"

#include "tinymat/dense/tools/dot.h"
#include "tinymat/dense/tools/norm2.h"

namespace tinymat
{
namespace DLA
{

#define SIGN( x, y ) (y) > 0 ? (x) : -(x)

//=====================================================================
template< int M, int N, class T >
void ApplyElementaryReflector( const VectorS< M, T >& v, const T& tau, MatrixS< M, N, T >& C )
/*
Based on the LAPACK routine DLARF

ApplyElementaryReflector applies a real elementary reflector H to a real m by n matrix
C, from either the left or the right. H is represented in the form

      H = I - tau * v * v**T

where tau is a real scalar and v is a real vector.

If tau = 0, then H is taken to be the unit matrix.
*/
{

  if ( tau == 0 )
    return;

  // Form  H * C

  MatrixS<1,N,Real> w = Transpose(v)*C;

  C -= tau*v*w;
}

//Compute sqrt( x'*x ) following the algorithm in LAPACK DNRM2
template< int M, class T >
T dnrm2(const VectorS< M, T >& x)
{
  T scale = 0;
  T ssq = 0;
  for ( int i = 0; i < M; i++ )
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


template< int M, class T >
void ElementaryReflector( T& alpha, VectorS< M, T >& x, T& tau )
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


#define INSTANCE_APPLY( M, N, T ) \
template void ApplyElementaryReflector<M,N,T>( const VectorS< M, T >& v, const T& tau, MatrixS< M, N, T >& C );

#define INSTANCE_REFLECTOR( M, T ) \
template void ElementaryReflector<M,T>( T& alpha, VectorS< M, T >& x, T& tau );


} //namespace DLA
} //namespace tinymat
