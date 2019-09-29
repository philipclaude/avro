// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXSYMS_DECOMPOSE_LDLT_H
#define MATRIXSYMS_DECOMPOSE_LDLT_H

#include "tools/SANSnumerics.h"     // Real
#include "MatrixS_Type.h"

namespace numpack 
{
namespace DLA
{

  //Perform the LDL^T decomposition of a symmetric matrix
  template< int M, class T >
  struct MatrixSymSLDLT
  {
    //cppcheck-suppress noExplicitConstructor
    MatrixSymSLDLT(const MatrixSymS<M, T>& MatrixExpr) : MatrixFac(MatrixExpr)
    {
      Decompose(MatrixFac);
    }

    static void Decompose( MatrixSymS<M, T>& Matrix );

    MatrixSymS<M, T> MatrixFac;
  };
#if 0
namespace Fixed
{

//Operator that allows for the syntax LDL = Decompose_LDLT(A)
  template< int M_, class T >
  struct MatrixDecompose_LDLT : public MatrixSType< MatrixDecompose_LDLT<M_,T>, true, false >
  {
    static const int M = M_;
    static const int N = M_;

    const MatrixSymS<M, T>& Matrix;

    explicit MatrixDecompose_LDLT( const MatrixSymS<M, T>& Matrix ) : Matrix(Matrix) {}

    inline void value(const Real sgn, MatrixS<M,N,T>& res) const
    {
      res = Matrix;
      MatrixSymSLDLT<M,T>::Decompose(res);
    }

    inline const MatrixDecompose_LDLT&
    operator+() const { return *this; }
  };

} //namespace Fixed

//=============================================================================
//Operator to generate a datatype to represent the LDLT decomposition of a symmetric matrix expression
  template< int M, class T >
  inline Fixed::MatrixDecompose_LDLT< M, T >
  Decompose_LDLT(const MatrixSymS<M, T>& Matrix)
  {
    return Fixed::MatrixDecompose_LDLT< M, T >(Matrix);
  }
#endif
  //A more traditional syntax Decompose_LDLT(A, LDL)
  template< int M, class T >
  inline void
  Decompose_LDLT(const MatrixSymS<M, T>& Matrix, MatrixSymS<M,T>& LD)
  {
    LD = Matrix;
    MatrixSymSLDLT< M, T >::Decompose( LD );
  }

} //namespace DLA
} //namespace numpack 


#endif //MATRIXSYMS_DECOMPOSE_LDLT_H
