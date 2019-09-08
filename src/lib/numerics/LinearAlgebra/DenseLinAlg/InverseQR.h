// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef INVERSEQR_H
#define INVERSEQR_H

#include "tools/SingularException.h"
#include "StaticSize/MatrixS_InverseQR.h"
#include "DynamicSize/MatrixD_InverseQR.h"

namespace SANS
{
namespace DLA
{

//=============================================================================
  //Operator to generate a datatype to represent the inverse (using QR without pivoting) of a matrix expression
  struct InverseQR
  {
    // Interfaces for computing a matrix inverse
    template< class Expr, bool useRF >
    static inline MatrixInverse< MatrixDQRSolver, Expr >
    Inverse(const MatrixDType<Expr, useRF>& MatrixExpr)
    {
      return MatrixInverse< MatrixDQRSolver, Expr >( MatrixExpr.cast() );
    }

    template< class Expr, bool useRF >
    static inline Fixed::MatrixInverse< MatrixSQRSolver, Expr >
    Inverse(const MatrixSType<Expr, useRF, true>& MatrixExpr)
    {
      return Fixed::MatrixInverse< MatrixSQRSolver, Expr >( MatrixExpr.cast() );
    }

    static inline Real
    Inverse(const Real& A)
    {
      SANS_ASSERT_NONSINGULAR(A);

      return Real(1)/A;
    }


    // Interface for performing a matrix solve
    template< class LExpr, bool useRFL, class RExpr, bool useRFR >
    static inline MatrixSolve< MatrixDQRSolver, typename LExpr::node_type, RExpr >
    Solve(const MatrixDType< LExpr, useRFL >& MatrixLExpr,
          const MatrixDType< RExpr, useRFR >& MatrixRExpr)
    {
      typedef typename LExpr::node_type T;
      return MatrixSolve< MatrixDQRSolver, T, RExpr >( MatrixLExpr, MatrixRExpr.cast());
    }

    template< class LExpr, bool useRFL, bool FullL, class RExpr, bool useRFR, bool FullR >
    static inline Fixed::MatrixSolve< MatrixSQRSolver, LExpr, RExpr >
    Solve(const MatrixSType<LExpr, useRFL, FullL>& MatrixLExpr,
          const MatrixSType<RExpr, useRFR, FullR>& MatrixRExpr)
    {
      return Fixed::MatrixSolve< MatrixSQRSolver, LExpr, RExpr >(MatrixLExpr.cast(), MatrixRExpr.cast());
    }

    static inline Real
    Solve(const Real& A, const Real& b)
    {
      SANS_ASSERT_NONSINGULAR(A);

      return b/A;
    }

    // Interface for generating a factorization class that can be used to solve a matrix multiple times
    template< class Expr, bool useRF >
    static inline MatrixFactorized< MatrixDQRSolver, typename Expr::node_type >
    Factorize(const MatrixDType< Expr, useRF >& MatrixExpr )
    {
      typedef typename Expr::node_type T;
      return MatrixFactorized< MatrixDQRSolver, T >( MatrixExpr );
    }

    template< class Expr, bool useRF, bool Full>
    static inline Fixed::MatrixFactorized< MatrixSQRSolver, Expr::M, Expr::N, typename Expr::Ttype >
    Factorize(const MatrixSType<Expr, useRF, Full>& MatrixExpr)
    {
      return Fixed::MatrixFactorized< MatrixSQRSolver, Expr::M, Expr::N, typename Expr::Ttype >(MatrixExpr);
    }

  };

} //namespace DLA
} //namespace SANS

#endif //INVERSEQR_H
