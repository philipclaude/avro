// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXD_QR_DECOMP_H
#define MATRIXD_QR_DECOMP_H

#include "tools/SANSnumerics.h"     // Real
#include "MatrixD_Type.h"
#include "VectorD.h"

namespace numpack 
{
namespace DLA
{

//Perform the LU decomposition with unit diagonal on U
template< class T >
struct MatrixDQR
{
  template<class Expr, bool useRF> //cppcheck-suppress noExplicitConstructor
  MatrixDQR(const MatrixDType<Expr, useRF>& MatrixExpr) : MatrixFac(MatrixExpr), tau(MatrixFac.n())
  {
    Decompose(MatrixFac, tau);
  }

  static void Decompose( MatrixDView< T >& Matrix, VectorDView<T>& tau );

  int m() const { return MatrixFac.n(); } //.n() is not a typo
  int n() const { return MatrixFac.n(); }
  int size() const { return m()*n(); }

  MatrixD<T> MatrixFac;
  VectorD<T> tau;
};
/*
//Operator that allows for the syntax QR = Decompose_QR(A)
template<class Expr>
struct MatrixDecompose_QR : public MatrixDType< MatrixDecompose_QR<Expr>, true >
{
  const Expr& MatrixExpr;
  MatrixDecompose_QR( const Expr& MatrixExpr ) : MatrixExpr(MatrixExpr) {}

  template< class T >
  inline void value(const Real sgn, MatrixDView<T>& res) const
  {
    res = MatrixExpr;
    MatrixDQR<T>::Decompose(res);
  }

  inline const MatrixDecompose_QR&
  operator+() const { return *this; }
  int m() const { return MatrixExpr.m(); }
  int n() const { return MatrixExpr.n(); }
  int size() const { return m()*n(); }
};

//=============================================================================
//Operator to generate a datatype to represent the QR decomposition without pivoting of a matrix expression
template< class Expr, bool useRF >
inline MatrixDecompose_QR< Expr >
Decompose_QR(const MatrixDType<Expr, useRF>& MatrixExpr)
{
  return MatrixDecompose_QR< Expr >(MatrixExpr.cast());
}
*/

//A more traditional syntax Decompose_QR(A, Q, R)
template< class Expr, bool useRF, class T >
inline void
Decompose_QR(const MatrixDType<Expr, useRF>& MatrixExpr, MatrixDView<T>& QR, VectorDView<T>& tau)
{
  QR = MatrixExpr;
  MatrixDQR< T >::Decompose( QR, tau );
}

} //namespace DLA
} //namespace numpack 


#endif //MATRIXD_QR_DECOMP_H
