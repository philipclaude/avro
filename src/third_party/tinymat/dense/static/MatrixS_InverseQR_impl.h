// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(MATRIXS_INVERSEQR_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "tinymat/dense/static/MatrixS.h"
#include "tinymat/dense/static/VectorS.h"
#include "tinymat/dense/InverseQR.h"

#include "ElementaryReflector.h"
#include "MatrixS_Decompose_QR.h"
#include "MatrixS_InverseQR.h"

//This computes computes a matrix inverse using QR decomposition without pivoting

//Based on LAPACK DGELS, DORM2R, and DTRSM

namespace tinymat
{
namespace DLA
{

template< int I, int NL >
struct QTrans
{
  template< int M, class T, class MatrixType >
  static void apply( MatrixS< M, NL, T >& A, const VectorS<NL,T>& tau, MatrixType& B )
  {
    static const int NR = MatrixType::N;

    //
    // B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS) (LAPACK DORM2R)
    //
    VectorS<M-I, T> v;
    v[0] = 1;
    for ( int i = 1; i < M-I; i++)
      v[i] = A(I+i,I);

    MatrixS<M-I,NR,T> C;
    for ( int i = 0; i < M-I; i++)
      for ( int j = 0; j < NR; j++)
        C(i,j) = B(I+i,j);

    //Apply H(i) to C(i:m,1:n) from the left
    ApplyElementaryReflector( v, tau( I ), C );

    for ( int i = 0; i < M-I; i++)
      for ( int j = 0; j < NR; j++)
        B(I+i,j) = C(i,j);

    QTrans<I+1,NL>::apply( A, tau, B );
  }
};

//Terminate the loop
template<int K>
struct QTrans<K,K>
{
  template< int M, int NL, class T, class MatrixType >
  static void apply( MatrixS< M, NL, T >& A, const VectorS<NL,T>& tau, MatrixType& B )
  {
  }
};

//-----------------------------------------------------------------------------
template< int M, int N, class T, class MatrixType >
template< int NR >
void MatrixSQRSolver<M, N, T, MatrixType>::Solve( const FactorType& Factorized, MatrixS< M, NR, T >& B, const Real& sgn, MatrixType& X )
{
  //BOOST_MPL_ASSERT_RELATION( M, >=, N );
  //BOOST_MPL_ASSERT_RELATION( N, ==, MatrixType::M );
  //BOOST_MPL_ASSERT_RELATION( NR, ==, MatrixType::N );

  //The algorithm from LAPACK requires temporary mods to AFac
  MatrixS< M, N, T >& AFac = const_cast<MatrixS< M, N, T >&>(Factorized.MatrixFac);
  const VectorS<N, T>& tau = Factorized.tau;

  // Least-Squares Problem min || A * X - B ||
  //
  // B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS) (DORM2R)
  //
  QTrans<0,N>::apply( AFac, tau, B );

  //
  // Form  B := alpha*inv( R )*B.
  //
  for (int j = N-1; j >= 0; --j)
  {
    for (int i = N-1; i > j; --i)
    {
      T factor = -AFac(j, i);
      B.axpy_rows(i, j, factor, 0, NR);
    }

    T invdiag = InverseQR::Inverse( AFac(j ,j) );
    B.scale_row(j, invdiag, 0, NR);
  }

  //Assign the solution to X
  for (int i = 0; i < N; i++)
    for (int j = 0; j < NR; j++)
      X(i,j) = sgn*B(i,j);

}


// Explicitly instantiate all datatypes used here
#define INSTANTIATE(M, N, NR, T, MatrixType ) \
template void MatrixSQRSolver< M, N, T, MatrixType >::Solve( const FactorType& Factorized, MatrixS< M, NR, T >& B, const Real& sgn, MatrixType& X );

#define MATRIXS( M, N, T ) MatrixS< M, N, T >

} //namespace DLA
} //namespace tinymat
