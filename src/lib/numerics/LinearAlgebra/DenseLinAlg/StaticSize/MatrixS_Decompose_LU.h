// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXS_DECOMPOSE_LU_H
#define MATRIXS_DECOMPOSE_LU_H

// Use boost static assert to show the integers in the compiler error messages.
// C++11 static_assert lacks this ability
//#include <boost/mpl/assert.hpp>

#include "tools/SANSnumerics.h"     // Real
#include "MatrixS_Type.h"
#include "LinearAlgebra/DenseLinAlg/tools/Identity.h"

namespace SANS
{
namespace DLA
{

//Perform the LU decomposition with unit diagonal on U
  template< int M, class T >
  struct MatrixSLU
  {
    template<class Expr, bool useRF> //cppcheck-suppress noExplicitConstructor
    MatrixSLU(const MatrixSType<Expr, useRF, true>& MatrixExpr) : MatrixFac(MatrixExpr)
    {
      Decompose(MatrixFac);
    }

    // special fast constructor for identity matrix ( i.e. no-op essentially)
    explicit MatrixSLU(const DLA::Identity& I) : MatrixFac(I) {}

    static void Decompose( MatrixS<M, M, T >& Matrix );

    MatrixS<M, M, T > MatrixFac;
  };

namespace Fixed
{

//Operator that allows for the syntax LU = Decompose_LU(A)
  template<class Expr>
  struct MatrixDecompose_LU : public MatrixSType< MatrixDecompose_LU<Expr>, true, true >
  {
    static const int M = Expr::M;
    static const int N = Expr::N;
    //BOOST_MPL_ASSERT_RELATION( Expr::M, ==, Expr::N );

    const Expr& MatrixExpr;

    explicit MatrixDecompose_LU( const Expr& MatrixExpr ) : MatrixExpr(MatrixExpr) {}

    template< class T >
    inline void value(const Real sgn, MatrixS<M,N,T>& res) const
    {
      res = MatrixExpr;
      MatrixSLU<M,T>::Decompose(res);
    }

    inline const MatrixDecompose_LU&
    operator+() const { return *this; }
  };

} //namespace Fixed

//=============================================================================
//Operator to generate a datatype to represent the LU decomposition without pivoting of a matrix expression
  template< class Expr, bool useRF >
  inline Fixed::MatrixDecompose_LU< Expr >
  Decompose_LU(const MatrixSType<Expr, useRF, true>& MatrixExpr)
  {
    return Fixed::MatrixDecompose_LU< Expr >(MatrixExpr.cast());
  }

//A more traditional syntax Decompose_LU(A, LU)
  template< int M, int N, class Expr, bool useRF, class T >
  inline void
  Decompose_LU(const MatrixSType<Expr, useRF, true>& MatrixExpr, MatrixS<M,N,T>& LU)
  {
    //BOOST_MPL_ASSERT_RELATION( Expr::M, ==, Expr::N );
    //BOOST_MPL_ASSERT_RELATION( M, ==, N );

    LU = MatrixExpr;
    MatrixSLU< M, T >::Decompose( LU );
  }

} //namespace DLA
} //namespace SANS


#endif //MATRIXS_DECOMPOSE_LU_H
