// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXD_LU_DECOMP_H
#define MATRIXD_LU_DECOMP_H

#include "tools/SANSnumerics.h"     // Real
#include "MatrixD_Type.h"

namespace SANS
{
namespace DLA
{

//Perform the LU decomposition with unit diagonal on U
template< class T >
struct MatrixDLU
{
  template<class Expr, bool useRF> //cppcheck-suppress noExplicitConstructor
  MatrixDLU(const MatrixDType<Expr, useRF>& MatrixExpr) : MatrixFac(MatrixExpr) { Decompose(MatrixFac); }

  static void Decompose( MatrixDView< T >& Matrix );

  int m() const { return MatrixFac.m(); }
  int n() const { return MatrixFac.n(); }
  typename MatrixD<T>::size_type size() const { return MatrixFac.size(); }

  MatrixD<T> MatrixFac;
};

//Operator that allows for the syntax LU = Decompose_LU(A)
template<class Expr>
struct MatrixDecompose_LU : public MatrixDType< MatrixDecompose_LU<Expr>, true >
{
  const Expr& MatrixExpr;
  explicit MatrixDecompose_LU( const Expr& MatrixExpr ) : MatrixExpr(MatrixExpr) {}

  template< class T >
  inline void value(const Real sgn, MatrixDView<T>& res) const
  {
    res = MatrixExpr;
    MatrixDLU<T>::Decompose(res);
  }

  inline const MatrixDecompose_LU&
  operator+() const { return *this; }
  int m() const { return MatrixExpr.m(); }
  int n() const { return MatrixExpr.n(); }
  int size() const { return m()*n(); }
};

//=============================================================================
//Operator to generate a datatype to represent the LU decomposition without pivoting of a matrix expression
template< class Expr, bool useRF >
inline MatrixDecompose_LU< Expr >
Decompose_LU(const MatrixDType<Expr, useRF>& MatrixExpr)
{
  return MatrixDecompose_LU< Expr >(MatrixExpr.cast());
}

//A more traditional syntax Decompose_LU(A, LU)
template< class Expr, bool useRF, class T >
inline void
Decompose_LU(const MatrixDType<Expr, useRF>& MatrixExpr, MatrixD<T>& LU)
{
  LU = MatrixExpr;
  MatrixDLU< T >::Decompose( LU );
}

} //namespace DLA
} //namespace SANS


#endif //MATRIXD_LU_DECOMP_H
