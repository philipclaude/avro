// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXD_ELEMENTARYREFLECTOR_H
#define MATRIXD_ELEMENTARYREFLECTOR_H

#include "MatrixD_Type.h"

namespace SANS
{
namespace DLA
{

template< class T, class T2 >
struct ElementaryReflector_impl
{
  //=============================================================================
  static void apply( const VectorDView< T >& v, const T& tau, MatrixDView< T2 >& C );
  /*
  Based on the LAPACK routine DLARF

  ApplyHouseholderMatrix applies a real elementary reflector H to a real m by n matrix
  C, from either the left or the right. H is represented in the form

        H = I - tau * v * v**T

  where tau is a real scalar and v is a real vector.

  If tau = 0, then H is taken to be the unit matrix.
  */

  //=============================================================================
  static void compute( T& alpha, VectorDView< T >& x, T& tau );
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
};

template<class T>
void ElementaryReflector( T& alpha, VectorDView< T >& x, T& tau )
{
  ElementaryReflector_impl<T,T>::compute(alpha, x, tau);
}

template<class T, class T2>
void ApplyElementaryReflector( const VectorDView< T >& v, const T& tau, MatrixDView< T2 >& C )
{
  ElementaryReflector_impl<T,T2>::apply(v, tau, C);
}

} //namespace DLA
} //namespace SANS

#endif //MATRIXD_ELEMENTARYREFLECTOR_H
