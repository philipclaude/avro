// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXBLOCK_2X2_DECOMPOSE_LU_H
#define MATRIXBLOCK_2X2_DECOMPOSE_LU_H

#include "tools/SANSnumerics.h"     // Real
#include "BlockLinAlg_Type.h"
#include "MatrixBlock_2x2.h"

namespace SANS
{
namespace BLA
{

//Perform the LU decomposition with unit diagonal on U
  template< class M00, class M01, class M10, class M11 >
  struct MatrixBlock_2x2_LU
  {
    template<class Expr> //cppcheck-suppress noExplicitConstructor
    MatrixBlock_2x2_LU(const BlockLinAlgType<Expr>& MatrixExpr) : MatrixFac(MatrixExpr.cast())
    {
      Decompose(MatrixFac);
    }

    static void Decompose( MatrixBlock_2x2<M00, M01, M10, M11>& Matrix );

    int m() const { return MatrixFac.m(); }
    int n() const { return MatrixFac.n(); }

    MatrixBlock_2x2<M00, M01, M10, M11> MatrixFac;
  };

  //Operator that allows for the syntax LU = Decompose_LU(A)
  template<class Expr>
  struct MatrixBlock_2x2_Decompose_LU : public BlockLinAlgType< MatrixBlock_2x2_Decompose_LU<Expr> >
  {
    explicit MatrixBlock_2x2_Decompose_LU( const Expr& MatrixExpr ) : MatrixExpr(MatrixExpr) {}

    template< class M00, class M01, class M10, class M11 >
    inline void value(const Real sgn, MatrixBlock_2x2<M00, M01, M10, M11>& res) const
    {
      res = MatrixExpr;
      MatrixBlock_2x2_LU<M00, M01, M10, M11>::Decompose(res);
    }

    inline const MatrixBlock_2x2_Decompose_LU&
    operator+() const { return *this; }

    int m() const { return MatrixExpr.m(); }
    int n() const { return MatrixExpr.n(); }
    int size() const { return m()*n(); }

    const Expr& MatrixExpr;
  };


//=============================================================================
//Operator to generate a datatype to represent the LU decomposition without pivoting of a matrix expression
  template< class Expr >
  inline MatrixBlock_2x2_Decompose_LU< Expr >
  Decompose_LU(const BlockLinAlgType<Expr>& MatrixExpr)
  {
    return MatrixBlock_2x2_Decompose_LU< Expr >(MatrixExpr.cast());
  }

  //A more traditional syntax Decompose_LU(A, LU)
  template< class Expr, class M00, class M01, class M10, class M11 >
  inline void
  Decompose_LU(const BlockLinAlgType<Expr>& MatrixExpr, MatrixBlock_2x2<M00, M01, M10, M11>& LU)
  {
    LU = MatrixExpr;
    MatrixBlock_2x2_LU<M00, M01, M10, M11>::Decompose( LU );
  }

//  //A more traditional syntax Decompose_LUP(A, LU)
//  //Note that LUP is only used for the inner block solves
//  template< class M00, class M01,
//            class M10, class M11, class Expr>
//  inline void
//  Decompose_LUP(const BlockLinAlgType<Expr>& MatrixExpr, MatrixBlock_2x2<M00, M01, M10, M11>& LU)
//  {
//  //    BOOST_MPL_ASSERT_RELATION( Expr::M, ==, Expr::N );
//  //    BOOST_MPL_ASSERT_RELATION( M, ==, N );
//
//      LU = MatrixExpr;
//      MatrixBlock_2x2_LU<M00, M01, M10, M11, DLA::InverseLUP>::Decompose( LU );
//  }

} //namespace BLA
} //namespace SANS


#endif //MATRIXBLOCK_2X2_DECOMPOSE_LU_H
