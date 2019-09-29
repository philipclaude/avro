// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXD_INVERSEQR_H
#define MATRIXD_INVERSEQR_H

#include "tools/SANSnumerics.h"     // Real
#include "MatrixD_Type.h"
#include "MatrixD_Inverse.h"
#include "MatrixD_Decompose_QR.h"

//Specialization of matrix inverse to solve a linear system using QR decomposition without pivoting

namespace numpack 
{
namespace DLA
{

// Structure for solving using QR without pivoting
  template< class T, class MatrixType >
  struct MatrixDQRSolver
  {
    //Convert the matrix to an QR decomposed matrix
    typedef MatrixDQR<T> FactorType;

    static void Solve(const FactorType& Factorized, MatrixDView< typename MatrixType::node_type >& B, const Real sgn, MatrixType& res );
  };

  //This series of static functions generate the appropriate temporary variables if necessary
  template< class T >
  struct InverseSolver< MatrixDQRSolver, T >
  {
    typedef typename MatrixDQRSolver<T, MatrixD<T> >::FactorType FactorType;

    // Inverse( A ) * (a + b)
    template< class RExpr, bool useRFR, class MatrixType >
    static void Solve( const FactorType& MatrixFac, const MatrixDType< RExpr, useRFR >& eR, const Real& sgn, MatrixType& res )
    {
      MatrixD< typename MatrixType::node_type > B(eR);
      Solve( MatrixFac, B, sgn, res );
    }

    // Inverse( A ) * b
    template< class MatrixType >
    static void Solve( const FactorType& MatrixFac, const MatrixDView< T >& MatrixR, const Real& sgn, MatrixType& res )
    {
      MatrixD< typename MatrixType::node_type > B(MatrixR);
      Solve(MatrixFac, B, sgn, res);
    }

    // Inverse( A )
    template< class MatrixType >
    static void Solve( const FactorType& MatrixFac, const Real& sgn, MatrixType& res )
    {
      MatrixD< typename MatrixType::node_type > B(MatrixFac.MatrixFac.m(), MatrixFac.MatrixFac.n());
      B = Identity();
      Solve(MatrixFac, B, sgn, res);
    }

    // This calls that actual inverse algorithm
    template< class MatrixType >
    static void Solve(const FactorType& MatrixFac, MatrixDView< typename MatrixType::node_type >& B, const Real& sgn, MatrixType& res )
    {
      MatrixDQRSolver<T, MatrixType>::Solve(MatrixFac, B, sgn, res);
    }
  };


} //namespace DLA
} //namespace numpack 

#endif //MATRIXD_INVERSEQR_H
