// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(MATRIXS_DECOMPOSE_QR_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "MatrixS_Decompose_QR.h"
#include "LinearAlgebra/DenseLinAlg/StaticSize/MatrixS.h"
#include "LinearAlgebra/DenseLinAlg/tools/dot.h"
#include "ElementaryReflector.h"

//Perform a QR decomposition without pivoting

namespace SANS
{
namespace DLA
{

//Based on LAPACK DGEQR2
//computes a QR factorization of a real M-by-N matrix A: A = Q * R

template<int I, int M, int N, class T>
void ComputeReflector(MatrixS< M, N, T >& A, VectorS<MIN(M,N),T>& tau)
{
  static const int xM = MAX(M-I-1,1);
  VectorS<xM, T> x;
  for ( int i = 0; i < xM; i++)
    x[i] = A(MIN(I+1+i,M-1),I);

  ElementaryReflector( A( I, I ), x, tau( I ) );

  for ( int i = 0; i < xM; i++)
    A(MIN(I+1+i,M-1),I) = x[i];
}


template<int I, int K, bool, bool>
struct Reflect
{
  template<int M, int N, class T>
  static void apply(MatrixS< M, N, T >& A, VectorS<MIN(M,N),T>& tau)
  {
    ComputeReflector<I>(A, tau);

    VectorS<M-I, T> v;
    v[0] = 1;
    for ( int i = 1; i < M-I; i++)
      v[i] = A(I+i,I);

    MatrixS<M-I,N-I-1,T> C;
    for ( int i = 0; i < M-I; i++)
      for ( int j = 0; j < N-I-1; j++)
        C(i,j) = A(I+i,I+j+1);

    ApplyElementaryReflector( v, tau( I ), C );

    for ( int i = 0; i < M-I; i++)
      for ( int j = 0; j < N-I-1; j++)
        A(I+i,I+j+1) = C(i,j);

    Reflect<I+1, K, I+1>=M-1, I+1>=N-1 >::apply(A,tau);
  }
};

template<int I, int K>
struct Reflect<I,K,true,false>
{
  template<int M, int N, class T>
  static void apply(MatrixS< M, N, T >& A, VectorS<MIN(M,N),T>& tau)
  {
    tau( I ) = 0;
    Reflect<I+1, K, I+1>=M-1, I+1>=N-1 >::apply(A,tau);
  }
};

template<int I, int K>
struct Reflect<I,K,false,true>
{
  template<int M, int N, class T>
  static void apply(MatrixS< M, N, T >& A, VectorS<MIN(M,N),T>& tau)
  {
    ComputeReflector<I>(A, tau);

    Reflect<I+1, K, I+1>=M-1, I+1>=N-1 >::apply(A,tau);
  }
};

//Terminate the loop
template<int K>
struct Reflect<K,K,false,true>
{
  template<int M, int N, class T>
  static void apply(MatrixS< M, N, T >& A, VectorS<MIN(M,N),T>& tau)
  {
    ComputeReflector<K>(A, tau);
  }
};

template<int K, bool b2>
struct Reflect<K,K,true,b2>
{
  template<int M, int N, class T>
  static void apply(MatrixS< M, N, T >& A, VectorS<MIN(M,N),T>& tau)
  {
    tau( K ) = 0;
  }
};

//-----------------------------------------------------------------------------
template< int M, int N, class T >
void MatrixSQR<M,N,T>::Decompose( MatrixS< M, N, T >& A, VectorS<MIN(M,N),T>& tau )
{
  //BOOST_MPL_ASSERT_RELATION( M, >=, N);

  const int K = MIN( M, N );

  Reflect<0, K-1, 0>=M-1, 0>=N-1 >::apply(A,tau);
  /*
  for ( int i = 0; i < K; i++ )
  {
    //Generate elementary reflector H(i) to annihilate A(i+1:m,i)
    if ( i >= m-1 )
      tau(i) = 0;
    else
    {
      VectorDView<T> x = A.subcol(std::min( i+1, M-1 ), i, M-i-1);
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
*/

}

} //namespace DLA
} //namespace SANS
