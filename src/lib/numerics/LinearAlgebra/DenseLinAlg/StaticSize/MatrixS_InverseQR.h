// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXS_INVERSEQR_H
#define MATRIXS_INVERSEQR_H

#include "tools/SANSnumerics.h"     // Real
#include "MatrixS_Type.h"
#include "MatrixS_Inverse.h"
#include "MatrixS_Decompose_QR.h"

//Specialization of matrix inverse to solve a linear system using QR decomposition without pivoting

namespace SANS
{
namespace DLA
{
  // Structure for solving using QR without pivoting
  template< int M, int N, class T, class MatrixType >
  struct MatrixSQRSolver
  {
    //Convert the matrix to an QR decomposed matrix
    typedef MatrixSQR<M, N, T> FactorType;

    template< int NR >
    static void Solve( const FactorType& Factorized, MatrixS< M, NR, T >& B, const Real& sgn, MatrixType& res );
  };

namespace Fixed
{
  //This series of static functions generate the appropriate temporary variables if necessary
  template< int M, int N, class T >
  struct InverseSolver< MatrixSQRSolver, M, N, T >
  {
    // Inverse( A ) * (a + b)
    template< class RExpr, bool useRFR, bool FullR, class MatrixType >
    static void Solve( const typename MatrixSQRSolver<M, N, T, MatrixType>::FactorType& MatrixFac,
                       const MatrixSType< RExpr, useRFR, FullR >& eR, const Real sgn, MatrixType& res )
    {
      MatrixS< RExpr::M, RExpr::N, typename MatrixType::Ttype > B(eR);
      Solve( MatrixFac, B, sgn, res );
    }

    // Inverse( A ) * b
    template< int NR, class MatrixType >
    static void Solve( const typename MatrixSQRSolver<M, N, T, MatrixType>::FactorType& MatrixFac,
                       const MatrixS< M, NR, T >& MatrixR, const Real sgn, MatrixType& res )
    {
      MatrixS< M, NR, T > B(MatrixR);
      Solve(MatrixFac, B, sgn, res);
    }

    // Inverse( A )
    template< class MatrixType >
    static void Solve( const typename MatrixSQRSolver<M, N, T, MatrixType>::FactorType& MatrixFac,
                       const Real& sgn, MatrixType& res )
    {
      MatrixS< MatrixType::M, MatrixType::N, typename MatrixType::Ttype > B;
      B = Identity();
      Solve(MatrixFac, B, sgn, res);
    }

    // This calls that actual inverse algorithm
    template< int NR, class MatrixType >
    static void Solve( const typename MatrixSQRSolver<M, N, T, MatrixType>::FactorType& MatrixFac,
                       MatrixS< M, NR, T >& B, const Real& sgn, MatrixType& res )
    {
      MatrixSQRSolver<M, N, T, MatrixType>::Solve(MatrixFac, B, sgn, res);
    }
  };

} //namespace Fixed
} //namespace DLA
} //namespace SANS

#endif //MATRIXS_INVERSEQR_H
