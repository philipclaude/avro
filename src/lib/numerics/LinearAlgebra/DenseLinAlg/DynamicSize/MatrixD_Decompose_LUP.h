// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXD_LUP_DECOMP_H
#define MATRIXD_LUP_DECOMP_H

#include "tools/SANSnumerics.h"     // Real
#include "MatrixD.h"
#include "VectorD.h"

namespace SANS
{
namespace DLA
{

//Perform the LU decomposition with unit diagonal on U
template< class T >
struct MatrixDLUP
{
  template<class Expr, bool useRF> //cppcheck-suppress noExplicitConstructor
  MatrixDLUP(const MatrixDType<Expr, useRF>& MatrixExpr) : MatrixFac(MatrixExpr), P(MatrixFac.m())
  {
    Decompose(MatrixFac, P);
  }

  static void Decompose( MatrixDView< T >& Matrix, VectorDView< int >& P );

  int m() const { return MatrixFac.m(); }
  int n() const { return MatrixFac.n(); }
  int size() const { return m()*n(); }

  MatrixD<T> MatrixFac;
  VectorD<int> P;
};

template<class Expr>
struct MatrixDecompose_LUP;

//This is a structure to store LU and P in the expression (LU, P) = Decompose_LUP(A);
//This structure represents the tuple (LU, P)
template< class T >
struct MatrixDecomposeTuple_LUP
{
  MatrixDView<T>& LU;
  VectorDView<int>& P;

  MatrixDecomposeTuple_LUP( MatrixDView<T>& LU, VectorDView<int>& P ) : LU(LU), P(P) {}

  template< class Expr >
  inline MatrixDecomposeTuple_LUP&
  operator=( const MatrixDecompose_LUP<Expr>& Decompose )
  {
    Decompose.value(*this);
    return *this;
  }
};

//This structure represents the function Decompose_LUP(A)
template<class Expr>
struct MatrixDecompose_LUP : public MatrixDType< MatrixDecompose_LUP<Expr>, true >
{
  typedef typename Expr::node_type node_type;
  typedef node_type T;

  const Expr& MatrixExpr;
  MatrixDecompose_LUP( const Expr& MatrixExpr ) : MatrixExpr(MatrixExpr) {}

  inline void value(MatrixDecomposeTuple_LUP<T>& res) const
  {
    res.LU = MatrixExpr;
    MatrixDLUP<T>::Decompose(res.LU, res.P);
  }

  inline const MatrixDecompose_LUP&
  operator+() const { return *this; }
  int m() const { return MatrixExpr.m(); }
  int n() const { return MatrixExpr.n(); }
  int size() const { return m()*n(); }
};


//=============================================================================
//Operator to generate a datatype to represent the LU decomposition with pivoting of a matrix expression
template< class Expr, bool useRF >
inline MatrixDecompose_LUP< Expr >
Decompose_LUP(const MatrixDType<Expr, useRF>& MatrixExpr)
{
  return MatrixDecompose_LUP< Expr >(MatrixExpr.cast());
}

//This operator generates the tuple (LU, P)
template<class T>
inline MatrixDecomposeTuple_LUP< T >
operator,( const MatrixDView<T>& LU, const VectorDView<int>& P)
{
  return MatrixDecomposeTuple_LUP<T>( const_cast<MatrixDView<T>&>(LU)
                                    , const_cast<VectorDView<int>&>(P) );
}


//=============================================================================
//A more classical way of calling the function Decompose_LUP(A, LU, P)
template< class Expr, bool useRF, class T >
inline void
Decompose_LUP(const MatrixDType<Expr, useRF>& Matrix, MatrixD<T>& LU, VectorD<int>& P)
{
  LU = Matrix;
  MatrixDLUP<T>::Decompose( LU, P );
}

} //namespace DLA
} //namespace SANS


#endif //MATRIXD_LUP_DECOMP_H
