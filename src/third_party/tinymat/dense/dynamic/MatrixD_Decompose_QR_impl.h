// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(MATRIXD_DECOMPOSE_QR_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "MatrixD_Decompose_QR.h"
#include "tinymat/dense/dynamic/MatrixD.h"
#include "tinymat/dense/static/MatrixS.h"
#include "tinymat/dense/tools/dot.h"
#include "ElementaryReflector.h"

//Perform a QR decomposition without pivoting

namespace tinymat 
{
namespace DLA
{

//Based on LAPACK DGEQR2
//computes a QR factorization of a real M-by-N matrix A: A = Q * R


//-----------------------------------------------------------------------------
template< class T >
void MatrixDQR<T>::Decompose( MatrixDView< T >& A, VectorDView< T >& tau )
{
  const int m = A.m();
  const int n = A.n();

  const int k = std::min( m, n );

  for ( int i = 0; i < k; i++ )
  {
    //Generate elementary reflector H(i) to annihilate A(i+1:m,i)
    if ( i >= m-1 )
      tau(i) = 0;
    else
    {
      VectorDView<T> x = A.subcol(std::min( i+1, m-1 ), i, m-i-1);
      ElementaryReflector( A( i, i ), x, tau( i ) );
    }
    if ( i < n-1 )
    {
      //Apply H(i) to A(i:m,i+1:n) from the left
      T aii = A( i, i );
      A( i, i ) = T(1);
      MatrixDView<T> C = A.sub(i,i+1,m-i,n-i-1);
      VectorDView<T> v = A.subcol(i, i, m-i);
      ApplyElementaryReflector( v, tau( i ), C );
      A( i, i ) = aii;
    }
  }


}

} //namespace DLA
} //namespace tinymat 
