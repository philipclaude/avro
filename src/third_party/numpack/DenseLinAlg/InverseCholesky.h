// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef INVERSECHOLESKY_H
#define INVERSECHOLESKY_H

#include "StaticSize/MatrixSymS_InverseCholesky.h"

namespace numpack 
{
namespace DLA
{

//=============================================================================
  //Operator to generate a datatype to represent the inverse (using LU without pivoting) of a matrix expression
  struct InverseCholesky
  {
    // Interfaces for computing a matrix inverse
    template< class Expr, bool useRF, bool Full >
    static inline Fixed::MatrixInverse< MatrixSymSCholeskySolver, Expr >
    Inverse(const MatrixSType<Expr, useRF, Full>& MatrixExpr)
    {
      return Fixed::MatrixInverse< MatrixSymSCholeskySolver, Expr >( MatrixExpr.cast() );
    }

    static inline Real
    Inverse(const Real& A)
    {
      return Real(1)/A;
    }

    // Interface for performing a matrix solve
    template< class LExpr, bool useRFL, bool FullL, class RExpr, bool useRFR, bool FullR >
    static inline Fixed::MatrixSolve< MatrixSymSCholeskySolver, LExpr, RExpr >
    Solve(const MatrixSType<LExpr, useRFL, FullL>& MatrixLExpr,
          const MatrixSType<RExpr, useRFR, FullR>& MatrixRExpr)
    {
      return Fixed::MatrixSolve< MatrixSymSCholeskySolver, LExpr, RExpr >(MatrixLExpr.cast(), MatrixRExpr.cast());
    }

    static inline Real
    Solve(const Real& A, const Real& b)
    {
      return b/A;
    }

    // Interface for generating a factorization class that can be used to solve a matrix multiple times
    template< class Expr, bool useRF, bool Full>
    static inline Fixed::MatrixFactorized< MatrixSymSCholeskySolver, Expr::M, Expr::N, typename Expr::Ttype >
    Factorize(const MatrixSType<Expr, useRF, Full>& MatrixExpr)
    {
      return Fixed::MatrixFactorized< MatrixSymSCholeskySolver, Expr::M, Expr::N, typename Expr::Ttype >(MatrixExpr);
    }

  };

} //namespace DLA
} //namespace numpack 

#endif //INVERSECHOLESKY_H
