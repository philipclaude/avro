// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef INVERSELDLT_H
#define INVERSELDLT_H

#include "tools/SingularException.h"

#include "static/MatrixSymS_InverseLDLT.h"

namespace tinymat 
{
namespace DLA
{

//=============================================================================
  //Operator to generate a datatype to represent the inverse (using LU without pivoting) of a matrix expression
  struct InverseLDLT
  {
    // Interfaces for computing a matrix inverse
    template< class Expr, bool useRF, bool Full >
    static inline Fixed::MatrixInverse< MatrixSymSLDLTSolver, Expr >
    Inverse(const MatrixSType<Expr, useRF, Full>& MatrixExpr)
    {
      return Fixed::MatrixInverse< MatrixSymSLDLTSolver, Expr >( MatrixExpr.cast() );
    }

    static inline Real
    Inverse(const Real& A)
    {
      SANS_ASSERT_NONSINGULAR(A);

      return Real(1)/A;
    }

    // Interface for performing a matrix solve
    template< class LExpr, bool useRFL, bool FullL, class RExpr, bool useRFR, bool FullR >
    static inline Fixed::MatrixSolve< MatrixSymSLDLTSolver, LExpr, RExpr >
    Solve(const MatrixSType<LExpr, useRFL, FullL>& MatrixLExpr,
          const MatrixSType<RExpr, useRFR, FullR>& MatrixRExpr)
    {
      return Fixed::MatrixSolve< MatrixSymSLDLTSolver, LExpr, RExpr >(MatrixLExpr.cast(), MatrixRExpr.cast());
    }

    static inline Real
    Solve(const Real& A, const Real& b)
    {
      SANS_ASSERT_NONSINGULAR(A);

      return b/A;
    }

    // Interface for generating a factorization class that can be used to solve a matrix multiple times
    template< class Expr, bool useRF, bool Full>
    static inline Fixed::MatrixFactorized< MatrixSymSLDLTSolver, Expr::M, Expr::N, typename Expr::Ttype >
    Factorize(const MatrixSType<Expr, useRF, Full>& MatrixExpr)
    {
      return Fixed::MatrixFactorized< MatrixSymSLDLTSolver, Expr::M, Expr::N, typename Expr::Ttype >(MatrixExpr);
    }

  };

} //namespace DLA
} //namespace tinymat 

#endif //INVERSELDLT_H
