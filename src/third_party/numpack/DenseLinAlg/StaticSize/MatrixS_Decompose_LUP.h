// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXS_LUP_DECOMP_H
#define MATRIXS_LUP_DECOMP_H

// Use boost static assert to show the integers in the compiler error messages.
// C++11 static_assert lacks this ability
//#include <boost/mpl/assert.hpp>

#include "tools/SANSnumerics.h"     // Real
#include "MatrixS_Type.h"
#include "VectorS.h"
#include "numpack/DenseLinAlg/tools/Identity.h"

namespace numpack 
{
namespace DLA
{

//Perform the LU decomposition with unit diagonal on U
  template< int M, class T >
  struct MatrixSLUP
  {
    template<class Expr, bool useRF> //cppcheck-suppress noExplicitConstructor
    MatrixSLUP(const MatrixSType<Expr, useRF, true>& MatrixExpr) : MatrixFac(MatrixExpr)
    {
     Decompose(MatrixFac, P);
    }

    // special fast constructor for identity matrix ( i.e. no-op essentially)
    explicit MatrixSLUP(const DLA::Identity& I);

    static void Decompose( MatrixS<M, M, T >& Matrix, VectorS<M, int>& P );

    MatrixS<M, M, T > MatrixFac;
    VectorS<M, int> P;
  };

namespace Fixed
{

  template<class Expr>
  struct MatrixDecompose_LUP;

  //This is a structure to store LU and P in the expression (LU, P) = Decompose_LUP(A);
  //This structure represents the tuple (LU, P)
  template< int M, class T >
  struct MatrixDecomposeTuple_LUP
  {
    MatrixS<M,M,T>& LU;
    VectorS<M,int>& P;

    MatrixDecomposeTuple_LUP( MatrixS<M,M,T>& LU, VectorS<M,int>& P ) : LU(LU), P(P) {}

    template< class Expr >
    inline MatrixDecomposeTuple_LUP&
    operator=( const MatrixDecompose_LUP<Expr>& Decompose )
    {
      Decompose.value(*this);
      return *this;
    }
  };

//Operator that allows for the syntax LU = Decompose_LU(A)
  template<class Expr>
  struct MatrixDecompose_LUP : public MatrixSType< MatrixDecompose_LUP<Expr>, true, true >
  {
    static const int M = Expr::M;
    static const int N = Expr::N;
    //BOOST_MPL_ASSERT_RELATION( Expr::M, ==, Expr::N );

    const Expr& MatrixExpr;

    explicit MatrixDecompose_LUP( const Expr& MatrixExpr ) : MatrixExpr(MatrixExpr) {}

    template< class T >
    inline void value(MatrixDecomposeTuple_LUP<M,T>& res) const
    {
      res.LU = MatrixExpr;
      MatrixSLUP<M,T>::Decompose(res.LU, res.P);
    }

    inline const MatrixDecompose_LUP&
    operator+() const { return *this; }
  };

} //namespace Fixed

//=============================================================================
//Operator to generate a datatype to represent the LU decomposition without pivoting of a matrix expression
  template< class Expr, bool useRF >
  inline Fixed::MatrixDecompose_LUP< Expr >
  Decompose_LUP(const MatrixSType<Expr, useRF, true>& MatrixExpr)
  {
    return Fixed::MatrixDecompose_LUP< Expr >(MatrixExpr.cast());
  }

  //This operator generates the tuple (LU, P)
  template<int M, class T>
  inline Fixed::MatrixDecomposeTuple_LUP< M, T >
  operator,( const MatrixS<M,M,T>& LU, const VectorS<M, int>& P)
  {
    return Fixed::MatrixDecomposeTuple_LUP<M, T>( const_cast<MatrixS<M,M,T>&>(LU),
                                                  const_cast<VectorS<M, int>&>(P) );
  }

//A more traditional syntax Decompose_LUP(A, LU, P)
  template< int M, int N, class Expr, bool useRF, class T >
  inline void
  Decompose_LUP(const MatrixSType<Expr, useRF, true>& MatrixExpr, MatrixS<M,N,T>& LU, VectorS<M, int>& P)
  {
    //BOOST_MPL_ASSERT_RELATION( Expr::M, ==, Expr::N );
    //BOOST_MPL_ASSERT_RELATION( M, ==, N );

    LU = MatrixExpr;
    MatrixSLUP< M, T >::Decompose( LU, P );
  }

} //namespace DLA
} //namespace numpack 


#endif //MATRIXS_LUP_DECOMP_H
