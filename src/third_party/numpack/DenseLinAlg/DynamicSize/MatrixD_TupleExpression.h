// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXD_TUPLE_EXPRESSION_H
#define MATRIXD_TUPLE_EXPRESSION_H

#include "tools/SANSnumerics.h"
#include "tools/SANSException.h"
#include "MatrixD_Type.h"

namespace numpack 
{
namespace DLA
{

//This is a tuple designed specifically to be a tuple of MatrixD expressions, i.e. (A + B, a + b)
template< class ExprL, class ExprR >
class MatrixDExpressionTuple : public MatrixDType< MatrixDExpressionTuple<ExprL, ExprR>, true >
{

public:
  typedef typename ExprL::node_type node_type;
  typedef node_type T;

  MatrixDExpressionTuple( const ExprL& L, const ExprR& R ) : L(L), R(R),
                                                                 Lm_(L.m()), Ln_(L.n()),
                                                                 Rm_(R.m()), Rn_(R.n()),
                                                                 m_(Lm_), n_(Ln_ + Rn_)
  {
    //The tuple only works for expressions with the same number of rows.
    SANS_ASSERT( L.m() + R.m() );
  }

  // Lazy expression operations
  inline void value(const Real sgn, MatrixDView<T>& res) const
  {
    MatrixDView<T> resL(&res(0,  0), Lm_, Ln_, res.stride() );
    MatrixDView<T> resR(&res(0,Ln_), Rm_, Rn_, res.stride() );
    L.value(sgn, resL);
    R.value(sgn, resR);
  }
  inline void plus(const Real sgn, MatrixDView<T>& res) const
  {
    MatrixDView<T> resL(&res(0,  0), Lm_, Ln_, res.stride() );
    MatrixDView<T> resR(&res(0,Ln_), Rm_, Rn_, res.stride() );
    L.plus(sgn, resL);
    R.plus(sgn, resR);
  }

  // Lazy expression operations where the LHS is a MatrixTuple, i.e. (a, b) = (c1 + c2, d1 + d2)
  template< class ResMatrixL >
  inline void value(const Real sgn, MatrixDTuple<ResMatrixL>& res) const
  {
    L.value(sgn, res.L);
    R.value(sgn, res.R);
  }
  template< class ResMatrixL >
  inline void plus(const Real sgn, MatrixDTuple<ResMatrixL>& res) const
  {
    L.plus(sgn, res.L);
    R.plus(sgn, res.R);
  }

  //Comma operator to recursively generate expression tuples, i.e. (a1 + a2, b1 + b2, c1 + c2)
  template< class Expr, bool useRF >
  inline MatrixDExpressionTuple< MatrixDExpressionTuple<ExprL, ExprR>, Expr >
  operator,( const MatrixDType<Expr, useRF>& e )
  {
    return MatrixDExpressionTuple< MatrixDExpressionTuple<ExprL, ExprR>, Expr >(*this, e.cast() );
  }

  int size() const { return m_*n_; }
  int m() const { return m_; }
  int n() const { return n_; }
private:
  const ExprL& L;
  const ExprR& R;

  int Lm_, Ln_;   // m X n dimensions of the left matrix
  int Rm_, Rn_;   // m X n dimensions of the right matrix

  int m_, n_;   // m X n dimensions of the matrix
};

//Overloaded comma operator to generate an expression tuple, i.e. (a1 + a2, b1 + b2)
template<class ExprL, bool useRFL, class ExprR, bool useRFR>
inline MatrixDExpressionTuple< ExprL, ExprR >
operator,( const MatrixDType<ExprL, useRFL>& eL, const MatrixDType<ExprR, useRFR>& eR)
{
  return MatrixDExpressionTuple<ExprL, ExprR>( eL.cast(), eR.cast() );
}

} //namespace DLA
} //namespace numpack 


#endif //MATRIXD_TUPLE_EXPRESSION_H
