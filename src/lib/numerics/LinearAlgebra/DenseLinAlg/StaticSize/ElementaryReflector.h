// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXS_ELEMENTARYREFLECTOR_H
#define MATRIXS_ELEMENTARYREFLECTOR_H

#include "MatrixS_Type.h"

namespace SANS
{
namespace DLA
{

//=============================================================================
template< int M, int N, class T >
void ApplyElementaryReflector( const VectorS< M, T >& v, const T& tau, MatrixS< M, N, T >& C );
/*
Based on the LAPACK routine DLARF

ApplyHouseholderMatrix applies a real elementary reflector H to a real m by n matrix
C, from either the left or the right. H is represented in the form

      H = I - tau * v * v**T

where tau is a real scalar and v is a real vector.

If tau = 0, then H is taken to be the unit matrix.
*/

//=============================================================================
template< int M, class T >
void ElementaryReflector( T& alpha, VectorS< M, T >& x, T& tau );
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

} //namespace DLA
} //namespace SANS

#endif //MATRIXS_ELEMENTARYREFLECTOR_H
