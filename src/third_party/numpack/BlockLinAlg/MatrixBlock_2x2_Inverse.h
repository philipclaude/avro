// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXBLOCK_2X2_INVERSE_H
#define MATRIXBLOCK_2X2_INVERSE_H

#include "tools/SANSnumerics.h"     // Real
#include "MatrixBlock_2x2.h"

#include "numpack/DenseLinAlg/tools/Identity.h"

//This is a general interfaces for adding matrix inverse operations. These classes are specialized for
//the specific type of inverse through the SolverType template argument

namespace numpack 
{
namespace BLA
{

  template< template<class, class, class, class, class> class SolverType, class M00, class M01, class M10, class M11 >
  struct InverseSolver
  {
    typedef MatrixBlock_2x2<M00, M01, M10, M11> MatrixType;
    typedef typename SolverType<M00, M01, M10, M11, MatrixType>::FactorType FactorType;

    // Inverse( A ) * b
    template< class ExpR, class MatrixType >
    static void Solve( const FactorType& MatrixFac,
                       const BlockLinAlgType<ExpR>& expr, const Real sgn, MatrixType& res )
    {
      res = expr;
      Solve(MatrixFac, sgn, res);
    }

    // This calls that actual inverse algorithm
    template< class MatrixType >
    static void Solve( const FactorType& MatrixFac,
                       const Real sgn, MatrixType& res )
    {
      SolverType<M00, M01, M10, M11, MatrixType>::Solve(MatrixFac, sgn, res);
    }
  };

  //Class for representing the inverse of a matrix
  template< template<class, class, class, class, class> class SolverType, class M00, class M01, class M10, class M11 >
  struct MatrixBlock_2x2_Inverse :
      public BlockLinAlgType< MatrixBlock_2x2_Inverse<SolverType, M00, M01, M10, M11> >
  {
    typedef MatrixBlock_2x2<M00, M01, M10, M11> MatrixType;
    typedef typename SolverType<M00, M01, M10, M11, MatrixType>::FactorType FactorType;

    const MatrixType& MatrixExpr; //The matrix expression that is solved

    // cppcheck-suppress noExplicitConstructor
    MatrixBlock_2x2_Inverse( const MatrixType& MatrixExpr ) : MatrixExpr(MatrixExpr) {}

    // A copy of MatrixExpr is created so that Inverse( A ) does not modify 'A'
    template<class RExpr>
    inline void value(const Real sgn, BlockLinAlgType<RExpr>& res) const
    {
      FactorType AFac(MatrixExpr);

      SANS_ASSERT(AFac.MatrixFac.m00.m() == res.cast().m00.m());
      SANS_ASSERT(AFac.MatrixFac.m00.n() == res.cast().m00.n());

      SANS_ASSERT(AFac.MatrixFac.m11.m() == res.cast().m11.m());
      SANS_ASSERT(AFac.MatrixFac.m11.n() == res.cast().m11.n());

      res.cast() = DLA::Identity();
      InverseSolver<SolverType, M00, M01, M10, M11>::Solve(AFac, sgn, res.cast());
    }

    inline const MatrixBlock_2x2_Inverse&
    operator+() const { return *this; }
  };


  //Represents the solution of a matrix inverse multiplied with a matrix expression, i.e. Solve(A,b);
  template< template<class, class, class, class, class> class SolverType,
            class M00, class M01, class M10, class M11, class RExpr>
  struct MatrixSolve : public BlockLinAlgType< MatrixSolve<SolverType, M00, M01, M10, M11, RExpr> >
  {
    typedef MatrixBlock_2x2<M00, M01, M10, M11> MatrixType;
    typedef typename SolverType<M00, M01, M10, M11, MatrixType>::FactorType FactorType;

    template<class Expr>
    MatrixSolve( const BlockLinAlgType<Expr>& Matrix, const RExpr& MatrixRExpr )
      : MatrixFac(Matrix), MatrixRExpr(MatrixRExpr) {}

    //This allows for the general expression x = s*Inverse( A + B )*(a + b)
    template<class Expr>
    inline void value(const Real sgn, BlockLinAlgType<Expr>& res) const
    {
      InverseSolver<SolverType, M00, M01, M10, M11>::Solve(MatrixFac, MatrixRExpr, sgn, res.cast());
    }

    inline const MatrixSolve&
    operator+() const { return *this; }

    int m() const { return MatrixFac.n(); }
    int n() const { return MatrixRExpr.n(); }
    int size() const { return m()*n(); }

  protected:
    const FactorType MatrixFac; //The factorized matrix that is solved
    const RExpr& MatrixRExpr;   //The matrix expression multiplying the inverse operation
  };

} //namespace BLA
} //namespace numpack 

#endif //MATRIXBLOCK_2X2_INVERSE_H
