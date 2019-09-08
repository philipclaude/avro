// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SRC_LINEARALGEBRA_BLOCKLINALG_BLOCKLINALG_MUL_H_
#define SRC_LINEARALGEBRA_BLOCKLINALG_BLOCKLINALG_MUL_H_

#include <type_traits>
//#include <boost/type_traits/is_arithmetic.hpp>

#include "tools/SANSnumerics.h"

#include "BlockLinAlg_Type.h"

namespace SANS
{
namespace BLA
{

// Multiplication involving BlockLinAlgType (vector, matrix or expression)

//=============================================================================
// Multiplication between a BlockLinAlg matrix or expression with a vector or matrix, including expressions
//
//-----------------------------------------------------------------------------
// This class represents multiplication between a matrix and a vector
// b = A*x
// no temporary variables are required
template<class Matrix_type, class Vector_type>
class OpMatMulVec : public BlockLinAlgType< OpMatMulVec<Matrix_type, Vector_type> >
{
public:
  OpMatMulVec(const Matrix_type& A, const Vector_type& x) :
    A_(A), x_(x)
  {
    SANS_ASSERT(A_.n()==x_.m());
    SANS_ASSERT(1==x_.n());
  }

  template<class Vector>
  inline void value(const Real sgn, BlockVectorType<Vector>& b) const
  {
    A_.mulVec_value(x_, sgn, b.cast());
  }
  template<class Vector>
  inline void plus(const Real sgn, BlockVectorType<Vector>& b) const
  {
    A_.mulVec_plus(x_, sgn, b.cast());
  }

  inline const OpMatMulVec& operator+() const { return *this; }

  int m() const { return A_.m(); }
  int n() const { return x_.n(); }

private:
  const Matrix_type& A_;
  const Vector_type& x_;
};

//-----------------------------------------------------------------------------
// This class represents multiplication between two 2x2 block matrices
// C = A*B
// no temporary variables are required
template<class MatrixL_type, class MatrixR_type>
class OpMatMul : public BlockLinAlgType< OpMatMul<MatrixL_type, MatrixR_type> >
{
public:
  OpMatMul(const MatrixL_type& A, const MatrixR_type& B) :
    A_(A), B_(B)
  {
    SANS_ASSERT(A_.n()==B_.m());
  }

  template<class Matrix>
  inline void value(const Real sgn, BlockMatrixType<Matrix>& C) const
  {
    A_.mulMat_value(B_, sgn, C.cast());
  }
  template<class Matrix>
  inline void plus(const Real sgn, BlockMatrixType<Matrix>& C) const
  {
    A_.mulMat_plus(B_, sgn, C.cast());
  }

  inline const OpMatMul& operator+() const { return *this; }

  int m() const { return A_.m(); }
  int n() const { return B_.n(); }

private:
  const MatrixL_type& A_;
  const MatrixR_type& B_;
};

//-----------------------------------------------------------------------------
// Overloaded operators for generating a representation of matrix multiplication (with vector or matrix) including expressions
template<class ExprL, class ExprR>
inline typename std::enable_if<isBlockMatrix<ExprL>::value && isBlockVector<ExprR>::value, OpMatMulVec<ExprL,ExprR>>::type
operator*(const BlockLinAlgType<ExprL>& eL, const BlockLinAlgType<ExprR>& eR)
{
  return OpMatMulVec<ExprL,ExprR>(eL.cast(), eR.cast());
}

template<class ExprL, class ExprR>
inline typename std::enable_if<isBlockMatrix<ExprL>::value && isBlockMatrix<ExprR>::value, OpMatMul<ExprL,ExprR>>::type
operator*(const BlockLinAlgType<ExprL>& eL, const BlockLinAlgType<ExprR>& eR)
{
  return OpMatMul<ExprL,ExprR>(eL.cast(), eR.cast());
}

//=============================================================================
// Multiplication between a BlockLinAlg expression and a scalar quantity
//
//-----------------------------------------------------------------------------
// This represents multiplication between a matrix or matrix expression and a scalar, i.e.
// M = s*(A + B)
// or
// M = (A + B)*s
// where s is a scalar quantity. The class also represents division with s by storing 1/s.
//
template<class Expr>
class OpMulScalar : public BlockLinAlgType< OpMulScalar<Expr> >
{
public:
  OpMulScalar(const Expr& e, const Real s) : e_(e), s_(s) {}

  template<class Vector>
  inline void value(const Real sgn, BlockVectorType<Vector>& res) const
  {
    e_.value(s_*sgn, res.cast());
  }

  template<class Vector>
  inline void plus(const Real sgn, BlockVectorType<Vector>& res) const
  {
    e_.plus(s_*sgn, res.cast());
  }

  inline const OpMulScalar& operator+() const { return *this; }

  int m() const { return e_.m(); }
  int n() const { return e_.n(); }

  const Expr& e_;
  const Real s_;
};

//-----------------------------------------------------------------------------
// Overloaded operators for generating a representation M of multiplication between a scalar and a matrix expression
// M = E*s
template<class Expr, typename T>
inline typename std::enable_if<std::is_arithmetic<T>::value, OpMulScalar<Expr>>::type // valid only if T is an arithmetic type
operator*(const BlockLinAlgType<Expr>& e, const T& s)
{
  return OpMulScalar<Expr>(e.cast(), s);
}
// M = E/s
template<class Expr, typename T>
inline typename std::enable_if<std::is_arithmetic<T>::value, OpMulScalar<Expr>>::type
operator/(const BlockLinAlgType<Expr>& e, const T& s)
{
  return OpMulScalar<Expr>(e.cast(), Real(1)/s);
}
// M = s*E
template<class Expr, typename T>
inline typename std::enable_if<std::is_arithmetic<T>::value, OpMulScalar<Expr>>::type
operator*(const T& s, const BlockLinAlgType<Expr>& e)
{
  return OpMulScalar<Expr>(e.cast(), s);
}
// M = -E
//Simple negation of an expression
template< class Expr>
inline const OpMulScalar<Expr>
operator-(BlockLinAlgType<Expr> const& e)
{
  return OpMulScalar<Expr>(e.cast(), -1);
}

//-----------------------------------------------------------------------------
// Overloaded operators a special case where matrix or matrix expression is multiplied from two sides, i.e. B = 2*A*3;
// This reduces the complexity of the expression tree and hence reduces code bloat
// M = (s1*E)*s2
template<class Expr, typename T>
inline typename std::enable_if<std::is_arithmetic<T>::value, OpMulScalar<Expr>>::type
operator*(const OpMulScalar<Expr>& MulScal, const T& s)
{
  return OpMulScalar<Expr>(MulScal.e_, MulScal.s_*s);
}
// M = (s1*E)/s2
template<class Expr, typename T>
inline typename std::enable_if< std::is_arithmetic<T>::value, OpMulScalar<Expr> >::type
operator/(const OpMulScalar<Expr>& MulScal, const T& s)
{
  return OpMulScalar<Expr>( MulScal.e_, MulScal.s_/s );
}
// M = s2*(s1*E)
template<class Expr, typename T>
inline typename std::enable_if< std::is_arithmetic<T>::value, OpMulScalar<Expr> >::type
operator*(const T& s, const OpMulScalar<Expr>& MulScal)
{
  return OpMulScalar<Expr>( MulScal.e_, MulScal.s_*s );
}
// M = -(s1*E)
template< class Expr >
inline const OpMulScalar<Expr>
operator-(OpMulScalar<Expr> const& Mul)
{
  return OpMulScalar<Expr>( Mul.e_, -1*Mul.s_ );
}

//=============================================================================
// BlockLinAlg matrix multiplication with vector or matrix that involves a scalar quantity
// such as s*A*x or A*x*s
//
template<class ExprL, class ExprR>
class OpMulFactor;

//-----------------------------------------------------------------------------
// This specialization translates (s*A)*x into s*(A*x)
// no temporary variables are required
template<class ExprL, class ExprR>
class OpMulFactor<OpMulScalar<ExprL>, ExprR> : public BlockLinAlgType< OpMulFactor<OpMulScalar<ExprL>, ExprR> >
{
public:

  OpMulFactor(const OpMulScalar<ExprL>& eL, const ExprR& eR) :
    eL_(eL), eR_(eR)
  {
    SANS_ASSERT(eL_.n() == eR_.m() );
  }

  template<class Vector>
  inline void value(const Real sgn, BlockVectorType<Vector>& res) const
  {
    res.cast() = eL_.s_*sgn*(eL_.e_*eR_);
  }

  template<class Vector>
  inline void plus(const Real sgn, BlockVectorType<Vector>& res) const
  {
    res.cast() += eL_.s_*sgn*(eL_.e_*eR_);
  }

  inline const OpMulFactor& operator+() const { return *this; }

  int m() const { return eL_.m(); }
  int n() const { return eR_.n(); }

private:
  const OpMulScalar<ExprL>& eL_;
  const ExprR& eR_;
};

//-----------------------------------------------------------------------------
// This specialization is translates A*(x*s) into (A*x)*s
// no temporary variables are required
template<class ExprL, class ExprR>
class OpMulFactor<ExprL, OpMulScalar<ExprR>> : public BlockLinAlgType< OpMulFactor<ExprL, OpMulScalar<ExprR>> >
{
public:

  OpMulFactor(const ExprL& eL, const OpMulScalar<ExprR>& eR) :
    eL_(eL), eR_(eR)
  {
    SANS_ASSERT(eL_.n() == eR_.m() );
  }

  template<class Vector>
  inline void value(const Real sgn, BlockVectorType<Vector>& res) const
  {
    res.cast() = eR_.s_*sgn*(eL_*eR_.e_);
  }

  template<class Vector>
  inline void plus(const Real sgn, BlockVectorType<Vector>& res) const
  {
    res.cast() += eR_.s_*sgn*(eL_*eR_.e_);
  }

  inline const OpMulFactor& operator+() const { return *this; }

  int m() const { return eL_.m(); }
  int n() const { return eR_.n(); }

private:
  const ExprL& eL_;
  const OpMulScalar<ExprR>& eR_;
};

//-----------------------------------------------------------------------------
// Overloaded operators for generating a representation of matrix multiplication including expressions
template<class ExprL, class ExprR>
inline OpMulFactor<OpMulScalar<ExprL>, ExprR>
operator*(const OpMulScalar<ExprL>& eL, const BlockLinAlgType<ExprR>& eR)
{
  return OpMulFactor<OpMulScalar<ExprL>, ExprR>(eL, eR.cast());
}

template<class ExprL, class ExprR>
inline OpMulFactor<ExprL, OpMulScalar<ExprR>>
operator*(const BlockLinAlgType<ExprL>& eL, const OpMulScalar<ExprR>& eR)
{
  return OpMulFactor<ExprL, OpMulScalar<ExprR>>(eL.cast(), eR);
}

} //namespace BLA
} //namespace SANS

#endif /* SRC_LINEARALGEBRA_BLOCKLINALG_BLOCKLINALG_MUL_H_ */
