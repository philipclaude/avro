// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(MATRIXD_INVERSELQR_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "numpack/DenseLinAlg/StaticSize/MatrixS.h"
#include "numpack/DenseLinAlg/StaticSize/VectorS.h"
#include "numpack/DenseLinAlg/DynamicSize/MatrixD.h"
#include "numpack/DenseLinAlg/DynamicSize/VectorD.h"
#include "numpack/DenseLinAlg/InverseQR.h"

#include "MatrixD_InverseQR.h"
#include "MatrixD_Decompose_QR.h"
#include "ElementaryReflector.h"

//This computes computes a matrix inverse using QR decomposition without pivoting

//Based on LAPACK DGELS, DORM2R, and DTRSM

namespace numpack 
{
namespace DLA
{

template< class T, class MatrixType >
void ApplyQTrans( MatrixDView< T >& A, const VectorD<T>& tau, MatrixType& B )
{
  const int m = A.m();
  const int k = A.n();
  const int n = B.n();

  //
  // B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS) (LAPACK DORM2R)
  //
  for ( int i = 0; i < k; i++ )
  {
    //Apply H(i) to C(i:m,1:n) from the left
    T aii = A( i, i );
    A( i, i ) = T(1);
    MatrixDView< typename MatrixType::node_type > C = B.sub(i, 0, m-i, n);
    VectorDView<T> v = A.subcol(i, i, m-i);
    ApplyElementaryReflector( v, tau( i ), C );
    A( i, i ) = aii;
  }
}


//-----------------------------------------------------------------------------
template< class T, class MatrixType >
void MatrixDQRSolver<T, MatrixType>::Solve(const FactorType& Factorized,
                                           MatrixDView< typename MatrixType::node_type >& B,
                                           const Real sgn, MatrixType& X )
{
  //The algorithm from LAPACK requires temporary mods to AFac
  MatrixD< T >& AFac = const_cast<MatrixD< T >&>(Factorized.MatrixFac);
  const VectorD<T>& tau = Factorized.tau;

  SANS_ASSERT( AFac.m() >= AFac.n() );
  SANS_ASSERT( AFac.m() == B.m() );
  SANS_ASSERT( AFac.n() == X.m() );

  const int n = AFac.n();
  const int nrhs = B.n();


  // Least-Squares Problem min || A * X - B ||
  //
  // B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS) (DORM2R)
  //
  ApplyQTrans( AFac, tau, B );

  //
  // Form  B := alpha*inv( R )*B.
  //
  for (int j = n-1; j >= 0; --j)
  {
    for (int i = n-1; i > j; --i)
    {
      T factor = -AFac(j, i);
      B.axpy_rows(i, j, factor, 0, nrhs);
    }

    T invdiag = InverseQR::Inverse( AFac(j ,j) );
    B.scale_row(j, invdiag, 0, nrhs);
  }

  X = sgn*B.sub(0, 0, n, nrhs);
}

} //namespace DLA
} //namespace numpack 
