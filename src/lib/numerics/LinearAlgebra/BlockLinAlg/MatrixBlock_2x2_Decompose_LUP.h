// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXBLOCK_2X2_DECOMPOSE_LUP_H
#define MATRIXBLOCK_2X2_DECOMPOSE_LUP_H

#include "tools/SANSnumerics.h"     // Real
#include "BlockLinAlg_Type.h"
#include "MatrixBlock_2x2.h"

namespace SANS
{
namespace BLA
{

//Perform the LU decomposition with unit diagonal on U, where LUP is used for the inner block solves
  template< class M00, class M01, class M10, class M11 >
  struct MatrixBlock_2x2_LUP
  {
    template<class Expr> //cppcheck-suppress noExplicitConstructor
    MatrixBlock_2x2_LUP(const BlockLinAlgType<Expr>& MatrixExpr) : MatrixFac(MatrixExpr.cast())
    {
      Decompose(MatrixFac);
    }

    static void Decompose( MatrixBlock_2x2<M00, M01, M10, M11>& Matrix );

    int m() const { return MatrixFac.m(); }
    int n() const { return MatrixFac.n(); }

    MatrixBlock_2x2<M00, M01, M10, M11> MatrixFac;
  };

  //Operator that allows for the syntax LU = Decompose_LUP(A)
  template<class Expr>
  struct MatrixBlock_2x2_Decompose_LUP : public BlockLinAlgType< MatrixBlock_2x2_Decompose_LUP<Expr> >
  {
    explicit MatrixBlock_2x2_Decompose_LUP( const Expr& MatrixExpr ) : MatrixExpr(MatrixExpr) {}

    template< class M00, class M01, class M10, class M11 >
    inline void value(const Real sgn, MatrixBlock_2x2<M00, M01, M10, M11>& res) const
    {
      res = MatrixExpr;
      MatrixBlock_2x2_LUP<M00, M01, M10, M11>::Decompose(res);
    }

    inline const MatrixBlock_2x2_Decompose_LUP&
    operator+() const { return *this; }

    int m() const { return MatrixExpr.m(); }
    int n() const { return MatrixExpr.n(); }
    int size() const { return m()*n(); }

    const Expr& MatrixExpr;
  };


//=============================================================================
//Operator to generate a datatype to represent the LU decomposition with pivoting of a matrix expression
  template< class Expr >
  inline MatrixBlock_2x2_Decompose_LUP< Expr >
  Decompose_LUP(const BlockLinAlgType<Expr>& MatrixExpr)
  {
    return MatrixBlock_2x2_Decompose_LUP< Expr >(MatrixExpr.cast());
  }

  //A more traditional syntax Decompose_LUP(A, LU)
  template< class Expr, class M00, class M01, class M10, class M11 >
  inline void
  Decompose_LUP(const BlockLinAlgType<Expr>& MatrixExpr, MatrixBlock_2x2<M00, M01, M10, M11>& LU)
  {
    LU = MatrixExpr;
    MatrixBlock_2x2_LUP<M00, M01, M10, M11>::Decompose( LU );
  }

} //namespace BLA
} //namespace SANS


#endif //MATRIXBLOCK_2X2_DECOMPOSE_LUP_H
