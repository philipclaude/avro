// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef INVERSELUP_H
#define INVERSELUP_H

#include "tools/SingularException.h"
#include "dynamic/MatrixD_InverseLUP.h"
#include "static/MatrixS_InverseLUP.h"
#include "numpack/block/MatrixBlock_2x2_InverseLUP.h"

namespace numpack 
{
namespace DLA
{

//=============================================================================
  //Operator to generate a datatype to represent the inverse (using LU with pivoting) of a matrix expression
  struct InverseLUP
  {
    // Interfaces for computing a matrix inverse
    template< class Expr, bool useRF >
    static inline MatrixInverse< MatrixDLUPSolver, Expr >
    Inverse(const MatrixDType<Expr, useRF>& MatrixExpr)
    {
      return MatrixInverse< MatrixDLUPSolver, Expr >( MatrixExpr.cast() );
    }

    template< class Expr, bool useRF >
    static inline Fixed::MatrixInverse< MatrixSLUPSolver, Expr >
    Inverse(const MatrixSType<Expr, useRF, true>& MatrixExpr)
    {
      return Fixed::MatrixInverse< MatrixSLUPSolver, Expr >( MatrixExpr.cast() );
    }

    template< class M00, class M01, class M10, class M11 >
    static inline BLA::MatrixBlock_2x2_Inverse< BLA::MatrixBlock_2x2_LUPSolver, M00, M01, M10, M11 >
    Inverse(const BLA::MatrixBlock_2x2<M00, M01, M10, M11>& Matrix)
    {
      return BLA::MatrixBlock_2x2_Inverse< BLA::MatrixBlock_2x2_LUPSolver, M00, M01, M10, M11 >( Matrix );
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
      SANS_ASSERT_NONSINGULAR(A);

      return Real(1)/A;
    }

    // Interface for performing a matrix solve
    template< class LExpr, bool useRFL, class RExpr, bool useRFR >
    static inline MatrixSolve< MatrixDLUPSolver, typename LExpr::node_type, RExpr >
    Solve(const MatrixDType< LExpr, useRFL >& MatrixLExpr,
          const MatrixDType< RExpr, useRFR >& MatrixRExpr)
    {
      typedef typename LExpr::node_type T;
      return MatrixSolve< MatrixDLUPSolver, T, RExpr >( MatrixLExpr, MatrixRExpr.cast() );
    }

    template< class LExpr, bool useRFL, bool FullL, class RExpr, bool useRFR, bool FullR >
    static inline Fixed::MatrixSolve< MatrixSLUPSolver, LExpr, RExpr >
    Solve(const MatrixSType<LExpr, useRFL, FullL>& MatrixLExpr,
          const MatrixSType<RExpr, useRFR, FullR>& MatrixRExpr)
    {
      return Fixed::MatrixSolve< MatrixSLUPSolver, LExpr, RExpr >( MatrixLExpr.cast(), MatrixRExpr.cast() );
    }

    template< class M00, class M01, class M10, class M11, class RExpr >
    static inline BLA::MatrixSolve< BLA::MatrixBlock_2x2_LUPSolver, M00, M01, M10, M11, RExpr >
    Solve(const BLA::MatrixBlock_2x2<M00, M01, M10, M11>& MatrixL,
          const BLA::blockType<RExpr>& MatrixRExpr)
    {
      return BLA::MatrixSolve< BLA::MatrixBlock_2x2_LUPSolver, M00, M01, M10, M11, RExpr >( MatrixL, MatrixRExpr.cast() );
    }

    static inline Real
    Solve(const Real& A, const Real& b)
    {
      SANS_ASSERT_NONSINGULAR(A);

      return b/A;
    }

    template<int N>
    static inline SurrealS<N,Real>
    Solve(const Real& A, const SurrealS<N,Real>& b)
    {
      SANS_ASSERT_NONSINGULAR(A);

      return b/A;
    }

    template<int N>
    static inline SurrealS<N,Real>
    Solve(const SurrealS<N,Real>& A, const SurrealS<N,Real>& b)
    {
      SANS_ASSERT_NONSINGULAR(A);

      return b/A;
    }

    // Interface for generating a factorization class that can be used to solve a matrix multiple times
    template< class Expr, bool useRF >
    static inline MatrixFactorized< MatrixDLUPSolver, typename Expr::node_type >
    Factorize(const MatrixDType< Expr, useRF >& MatrixExpr )
    {
      typedef typename Expr::node_type T;
      return MatrixFactorized< MatrixDLUPSolver, T >( MatrixExpr );
    }

    template< class Expr, bool useRF, bool Full>
    static inline Fixed::MatrixFactorized< MatrixSLUPSolver, Expr::M, Expr::N, typename Expr::Ttype >
    Factorize(const MatrixSType<Expr, useRF, Full>& MatrixExpr)
    {
      return Fixed::MatrixFactorized< MatrixSLUPSolver, Expr::M, Expr::N, typename Expr::Ttype >(MatrixExpr);
    }

  };

} //namespace DLA
} //namespace numpack 

#endif //INVERSE_H
