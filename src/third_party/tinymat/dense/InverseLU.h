// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef INVERSELU_H
#define INVERSELU_H

#include "tools/SingularException.h"
#include "dynamic/MatrixD_InverseLU.h"
#include "static/MatrixS_InverseLU.h"

namespace tinymat
{
namespace DLA
{

//=============================================================================
  //Operator to generate a datatype to represent the inverse (using LU without pivoting) of a matrix expression
  struct InverseLU
  {
    // Interfaces for computing a matrix inverse
    template< class Expr, bool useRF >
    static inline MatrixInverse< MatrixDLUSolver, Expr >
    Inverse(const MatrixDType<Expr, useRF>& MatrixExpr)
    {
      return MatrixInverse< MatrixDLUSolver, Expr >( MatrixExpr.cast() );
    }

    template< class Expr, bool useRF >
    static inline Fixed::MatrixInverse< MatrixSLUSolver, Expr >
    Inverse(const MatrixSType<Expr, useRF, true>& MatrixExpr)
    {
      return Fixed::MatrixInverse< MatrixSLUSolver, Expr >( MatrixExpr.cast() );
    }

    static inline Real
    Inverse(const Real& A)
    {
      SANS_ASSERT_NONSINGULAR(A);

      return Real(1)/A;
    }

    template<int N>
    static inline SurrealS<N,Real>
    Inverse(const SurrealS<N,Real>& A)
    {
      SANS_ASSERT_NONSINGULAR(A.value());

      return Real(1)/A;
    }

    // Interface for performing a matrix solve
    template< class LExpr, bool useRFL, class RExpr, bool useRFR >
    static inline MatrixSolve< MatrixDLUSolver, typename LExpr::node_type, RExpr >
    Solve(const MatrixDType< LExpr, useRFL >& MatrixLExpr,
          const MatrixDType< RExpr, useRFR >& MatrixRExpr)
    {
      typedef typename LExpr::node_type T;
      return MatrixSolve< MatrixDLUSolver, T, RExpr >( MatrixLExpr, MatrixRExpr.cast() );
    }

    template< class LExpr, bool useRFL, bool FullL, class RExpr, bool useRFR, bool FullR >
    static inline Fixed::MatrixSolve< MatrixSLUSolver, LExpr, RExpr >
    Solve(const MatrixSType<LExpr, useRFL, FullL>& MatrixLExpr,
          const MatrixSType<RExpr, useRFR, FullR>& MatrixRExpr)
    {
      return Fixed::MatrixSolve< MatrixSLUSolver, LExpr, RExpr >(MatrixLExpr.cast(), MatrixRExpr.cast());
    }

    static inline Real
    Solve(const Real& A, const Real& b)
    {
      SANS_ASSERT_NONSINGULAR(A);

      return b/A;
    }

    template<int N>
    static inline SurrealS<N,Real>
    Solve(const SurrealS<N,Real>& A, const SurrealS<N,Real>& b)
    {
      SANS_ASSERT_NONSINGULAR(A.value());

      return b/A;
    }

    // Interface for generating a factorization class that can be used to solve a matrix multiple times
    template< class Expr, bool useRF >
    static inline MatrixFactorized< MatrixDLUSolver, typename Expr::node_type >
    Factorize(const MatrixDType< Expr, useRF >& MatrixExpr )
    {
      typedef typename Expr::node_type T;
      return MatrixFactorized< MatrixDLUSolver, T >( MatrixExpr );
    }

    template< class Expr, bool useRF, bool Full>
    static inline Fixed::MatrixFactorized< MatrixSLUSolver, Expr::M, Expr::N, typename Expr::Ttype >
    Factorize(const MatrixSType<Expr, useRF, Full>& MatrixExpr)
    {
      return Fixed::MatrixFactorized< MatrixSLUSolver, Expr::M, Expr::N, typename Expr::Ttype >(MatrixExpr);
    }

  };

} //namespace DLA
} //namespace tinymat

#endif //INVERSELU_H
