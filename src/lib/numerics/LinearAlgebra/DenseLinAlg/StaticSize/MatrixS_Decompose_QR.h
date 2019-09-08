// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXS_QR_DECOMP_H
#define MATRIXS_QR_DECOMP_H

#include "tools/SANSnumerics.h"     // Real
#include "tools/minmax.h"
#include "MatrixS_Type.h"
#include "VectorS.h"

namespace SANS
{
namespace DLA
{

//Perform the QR decomposition
template< int M, int N, class T >
struct MatrixSQR
{
  template<class Expr, bool useRF> //cppcheck-suppress noExplicitConstructor
  MatrixSQR(const MatrixSType<Expr, useRF, true>& MatrixExpr) : MatrixFac(MatrixExpr)
  {
    Decompose(MatrixFac, tau);
  }

  static void Decompose( MatrixS< M, N, T >& A, VectorS<MIN(M,N),T>& tau );

  MatrixS<M, N, T > MatrixFac;
  VectorS<MIN(M,N),T> tau;
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

//A more traditional syntax Decompose_QR(A, Q, tau)
template< class Expr, bool useRF, bool MatrixFull, int M, int N, class T >
inline void
Decompose_QR(const MatrixSType<Expr, useRF, MatrixFull>& MatrixExpr, MatrixS<M,N,T>& QR, VectorS<MIN(M,N),T>& tau)
{
  QR = MatrixExpr;
  MatrixSQR< M, N, T >::Decompose( QR, tau );
}

} //namespace DLA
} //namespace SANS


#endif //MATRIXS_QR_DECOMP_H
