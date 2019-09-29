// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXSYMS_DECOMPOSE_CHOLESKY_H
#define MATRIXSYMS_DECOMPOSE_CHOLESKY_H

#include "tools/SANSnumerics.h"     // Real
#include "MatrixS_Type.h"

namespace numpack 
{
namespace DLA
{

//Perform the Cholesky decomposition of a symmetric positive definite matrix
  template< int M, class T >
  struct MatrixSymSCholesky
  {
    //cppcheck-suppress noExplicitConstructor
    MatrixSymSCholesky(const MatrixSymS<M, T>& MatrixExpr) : MatrixFac(MatrixExpr)
    {
      Decompose(MatrixFac);
    }

    static void Decompose( MatrixSymS<M, T>& Matrix );

    MatrixSymS<M, T> MatrixFac;
  };
#if 0
namespace Fixed
{

//Operator that allows for the syntax LL = Decompose_Cholesky(A)
  template< int M_, class T >
  struct MatrixDecompose_Cholesky : public MatrixSType< MatrixDecompose_Cholesky<M_,T>, true, false >
  {
    static const int M = M_;
    static const int N = M_;

    const MatrixSymS<M, T>& Matrix;

    explicit MatrixDecompose_Cholesky( const MatrixSymS<M, T>& Matrix ) : Matrix(Matrix) {}

    inline void value(const Real sgn, MatrixS<M,N,T>& res) const
    {
      res = Matrix;
      MatrixSymSCholesky<M,T>::Decompose(res);
    }

    inline const MatrixDecompose_Cholesky&
    operator+() const { return *this; }
  };

} //namespace Fixed

//=============================================================================
//Operator to generate a datatype to represent the Cholesky decomposition of a symmetric matrix expression
  template< int M, class T >
  inline Fixed::MatrixDecompose_Cholesky< M, T >
  Decompose_Cholesky(const MatrixSymS<M, T>& Matrix)
  {
    return Fixed::MatrixDecompose_Cholesky< M, T >(Matrix);
  }
#endif
//A more traditional syntax Decompose_Cholesky(A, LL)
  template< int M, class T >
  inline void
  Decompose_Cholesky(const MatrixSymS<M, T>& Matrix, MatrixSymS<M,T>& LL)
  {
    LL = Matrix;
    MatrixSymSCholesky< M, T >::Decompose( LL );
  }

} //namespace DLA
} //namespace numpack 


#endif //MATRIXSYMS_DECOMPOSE_CHOLESKY_H
