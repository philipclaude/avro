// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SPARSELINALG_MUL_H
#define SPARSELINALG_MUL_H

#include <boost/type_traits/is_arithmetic.hpp>
#include <memory> //shared_ptr

#include "SparseLinAlg_Type.h"

#include "SparseMatrix_CRS.h"

namespace SANS
{
namespace SLA
{

template<class , class >
class OpMul;

//-----------------------------------------------------------------------------
// This class represents multiplication between a matrix and a vector
// b = A*x
// no temporary variables are required
template<class Matrix_type, class TV>
class OpMatMulVec : public SparseLinAlgType< OpMatMulVec< Matrix_type, TV >, true >
{
public:
  typedef typename Matrix_type::Ttype Ttype;

  OpMatMulVec(const Matrix_type& A, const SparseVector<TV>& x) : A(A), x(x) {}

  template<class T>
  inline void value(const Real sgn, SparseVector<T>& b) const
  {
    A.mulVec_value( x, sgn ,b );
  }
  template<class T>
  inline void plus(const Real sgn, SparseVector<T>& b) const
  {
    A.mulVec_plus( x, sgn ,b );
  }

  inline const OpMatMulVec&
  operator+() const { return *this; }
  int m() const { return A.m(); }
private:
  const Matrix_type& A;
  const SparseVector<TV>& x;
};

//This is a specialization of the general OpMul for multiplication between CRS matrix and a vector
template<class TM, class TV>
class OpMul< SparseMatrix_CRS<TM>, SparseVector<TV> > : public OpMatMulVec< SparseMatrix_CRS<TM>, TV >
{
public:
  typedef OpMatMulVec< SparseMatrix_CRS<TM>, TV > Base_type;
  OpMul(const SparseMatrix_CRS<TM>& A, const SparseVector<TV>& x) : Base_type::OpMatMulVec(A, x) {}
};


//-----------------------------------------------------------------------------
// Operator for generating an OpMul representation of matrix multiplication including expressions
template<class ExprL, bool useRFL, class ExprR, bool useRFR>
inline OpMul<ExprL, ExprR>
operator*(const SparseLinAlgType<ExprL, useRFL>& eL, const SparseLinAlgType<ExprR, useRFR>& eR)
{
  return OpMul<ExprL, ExprR>( eL.cast(), eR.cast() );
}

template<class TM, class ExprR, bool useRFR>
inline OpMul< SparseMatrix_CRS<TM>, ExprR >
operator*(const std::shared_ptr< SparseMatrix_CRS<TM> >& eL, const SparseLinAlgType<ExprR, useRFR>& eR)
{
  return OpMul< SparseMatrix_CRS<TM>, ExprR>( *eL.get(), eR.cast() );
}

// Multiplication between a SparseLinAlg expression and a scalar quantity

//-----------------------------------------------------------------------------
// This represents multiplication between a matrix or matrix expression and a scalar, i.e.
// M = s*(A + B)
// or
// M = (A + B)*s
// where s is a scalar quantity. The class also represents division with s by storing 1/s.
template<class Expr, bool useRF>
class OpMulScalar;

template<class Expr>
class OpMulScalar<Expr, true> : public SparseLinAlgType< OpMulScalar<Expr, true>, true >
{
public:
  typedef typename Expr::Ttype Ttype;

  OpMulScalar(const Expr& e, const Real s) : e(e), s(s) {}

  template<class TV>
  inline void value(const Real sgn, SparseVector<TV>& res) const
  {
    e.value(s*sgn, res);
  }
  template<class TV>
  inline void plus(const Real sgn, SparseVector<TV>& res) const
  {
    e.plus(s*sgn, res);
  }

  inline const OpMulScalar&
  operator+() const { return *this; }
  int m() const { return e.m(); }

  const Expr& e;
  const Real s;
};

template<class Expr>
class OpMulScalar<Expr, false> : public SparseLinAlgType< OpMulScalar<Expr, false>, false >
{
public:
  typedef typename Expr::Ttype Ttype;

  OpMulScalar(const Expr& e, const Real s) : e(e), s(s) {}

  //Element-wise expression
  inline Ttype operator[](const int& i) const { return s*e[i]; }

  inline const OpMulScalar&
  operator+() const { return *this; }
  int m() const { return e.m(); }

  const Expr& e;
  const Real s;
};


//=============================================================================
// Overloaded operators to represent multiplication between a scalar and a matrix expression
template<class Expr, bool useRF, typename T>
inline typename std::enable_if< boost::is_arithmetic<T>::value, OpMulScalar<Expr, useRF> >::type
operator*(const SparseLinAlgType<Expr, useRF>& e, const T& s)
{
  return OpMulScalar<Expr, useRF>( e.cast(), s );
}

template<class Expr, bool useRF, typename T>
inline typename std::enable_if< boost::is_arithmetic<T>::value, OpMulScalar<Expr, useRF> >::type
operator/(const SparseLinAlgType<Expr, useRF>& e, const T& s)
{
  return OpMulScalar<Expr, useRF>( e.cast(), Real(1)/s );
}

template<class Expr, bool useRF, typename T>
inline typename std::enable_if< boost::is_arithmetic<T>::value, OpMulScalar<Expr, useRF> >::type
operator*(const T& s, const SparseLinAlgType<Expr, useRF>& e)
{
  return OpMulScalar<Expr, useRF>( e.cast(), s );
}

//Simple negation of an expression
template< class Expr, bool useRF >
inline const OpMulScalar<Expr, useRF>
operator-(SparseLinAlgType<Expr, useRF> const& e)
{
  return OpMulScalar<Expr, useRF>( e.cast(), -1 );
}

//=============================================================================
// This is a special case when multiplies scalars are multiplying from two sides, i.e. B = 2*A*2;
// This reduces the complexity of the expression tree and hence reduces code bloat
template<class Expr, bool useRF, typename T>
inline typename std::enable_if< boost::is_arithmetic<T>::value, OpMulScalar<Expr, useRF> >::type
operator*(const OpMulScalar<Expr, useRF>& MulScal, const T& s)
{
  return OpMulScalar<Expr, useRF>( MulScal.e, MulScal.s*s );
}

template<class Expr, bool useRF, typename T>
inline typename std::enable_if< boost::is_arithmetic<T>::value, OpMulScalar<Expr, useRF> >::type
operator/(const OpMulScalar<Expr, useRF>& MulScal, const T& s)
{
  return OpMulScalar<Expr, useRF>( MulScal.e, MulScal.s/s );
}

template<class Expr, bool useRF, typename T>
inline typename std::enable_if< boost::is_arithmetic<T>::value, OpMulScalar<Expr, useRF> >::type
operator*(const T& s, const OpMulScalar<Expr, useRF>& MulScal)
{
  return OpMulScalar<Expr, useRF>( MulScal.e, MulScal.s*s );
}

template< class Expr, bool useRF >
inline const OpMulScalar<Expr, useRF>
operator-(OpMulScalar<Expr, useRF> const& Mul)
{
  return OpMulScalar<Expr, useRF>( Mul.e, -1*Mul.s );
}



template<class ExprL, class ExprR>
class OpMulFactor;

//-----------------------------------------------------------------------------
// This specialization is for multiplication that involve a scalar and translates
// (s*A)*x into s*(A*x)
// no temporary variables are required
template<class ExprL, bool useRFL, class ExprR>
class OpMulFactor< OpMulScalar<ExprL, useRFL>, ExprR > : public SparseLinAlgType< OpMulFactor< OpMulScalar<ExprL, useRFL>, ExprR >, true >
{
public:
  typedef typename ExprL::Ttype Ttype;

  OpMulFactor(const OpMulScalar<ExprL, useRFL>& eL, const ExprR& eR) : eL(eL), eR(eR) { SANS_ASSERT(eL.m() == eR.m() ); }

  template<class TV>
  inline void value(const Real sgn, SparseVector<TV>& res) const
  {
    res = eL.s*sgn*(eL.e*eR);
  }
  template<class TV>
  inline void plus(const Real sgn, SparseVector<TV>& res) const
  {
    res += eL.s*sgn*(eL.e*eR);
  }

  inline const OpMulFactor&
  operator+() const { return *this; }
  int m() const { return eR.m(); }
private:
  const OpMulScalar<ExprL, useRFL>& eL;
  const ExprR& eR;
};

//-----------------------------------------------------------------------------
// This specialization is for multiplication that involve a scalar and translates
// A*(x*s) into (A*x)*s
// no temporary variables are required
template<class ExprL, class ExprR, bool useRFR>
class OpMulFactor< ExprL, OpMulScalar<ExprR, useRFR> > : public SparseLinAlgType< OpMulFactor< ExprL, OpMulScalar<ExprR, useRFR> >, true >
{
public:
  typedef typename ExprL::Ttype Ttype;

  OpMulFactor(const ExprL& eL, const OpMulScalar<ExprR, useRFR>& eR) : eL(eL), eR(eR) { SANS_ASSERT(eL.m() == eR.m() ); }

  template<class TV>
  inline void value(const Real sgn, SparseVector<TV>& res) const
  {
    res = eR.s*sgn*(eL*eR.e);
  }
  template<class TV>
  inline void plus(const Real sgn, SparseVector<TV>& res) const
  {
    res += eR.s*sgn*(eL*eR.e);
  }

  inline const OpMulFactor&
  operator+() const { return *this; }
  int m() const { return eR.m(); }
private:
  const ExprL& eL;
  const OpMulScalar<ExprR, useRFR>& eR;
};


//-----------------------------------------------------------------------------
// Operator for generating an OpMul representation of matrix multiplication including expressions
template<class ExprL, bool useRFL, class ExprR, bool useRFR>
inline OpMulFactor< OpMulScalar<ExprL, useRFL>, ExprR>
operator*(const OpMulScalar<ExprL, useRFL>& eL, const SparseLinAlgType<ExprR, useRFR>& eR)
{
  return OpMulFactor< OpMulScalar<ExprL, useRFL>, ExprR>( eL, eR.cast() );
}

template<class ExprL, bool useRFL, class ExprR, bool useRFR>
inline OpMulFactor<ExprL, OpMulScalar<ExprR, useRFR> >
operator*(const SparseLinAlgType<ExprL, useRFL>& eL, const OpMulScalar<ExprR, useRFR>& eR)
{
  return OpMulFactor<ExprL, OpMulScalar<ExprR, useRFR> >( eL.cast(), eR );
}


//-----------------------------------------------------------------------------
// Dummy class for NonZeroPattern multiplication
template<class TM, class TV>
class OpMatMulVec<SparseNonZeroPattern<TM>, TV> : public SparseLinAlgType< OpMatMulVec< SparseNonZeroPattern<TM>, TV >, true >
{
public:
  typedef SparseNonZeroPattern<TM> Matrix_type;
  typedef TM Ttype;

  OpMatMulVec(const Matrix_type& A, const SparseVector<TV>& x) : A(A) {}

  template<class T>
  inline void value(const Real sgn, SparseVector<T>& b) const
  {
  }
  template<class T>
  inline void plus(const Real sgn, SparseVector<T>& b) const
  {
  }

  inline const OpMatMulVec&
  operator+() const { return *this; }
  int m() const { return A.m(); }
private:
  const Matrix_type& A;
};

//This is a specialization of the general OpMul for multiplication between CRS matrix and a vector
template<class TM, class TV>
class OpMul< SparseNonZeroPattern<TM>, SparseVector<TV> > : public OpMatMulVec< SparseNonZeroPattern<TM>, TV >
{
public:
  typedef OpMatMulVec< SparseNonZeroPattern<TM>, TV > Base_type;
  OpMul(const SparseNonZeroPattern<TM>& A, const SparseVector<TV>& x) : Base_type::OpMatMulVec(A, x) {}
};


//-----------------------------------------------------------------------------
// Operator for generating an OpMul representation of matrix multiplication including expressions
template<class TM, class TV>
inline OpMul<SparseNonZeroPattern<TM>, SparseVector<TV>>
operator*(const SparseNonZeroPattern<TM>& A, const SparseVector<TV>& x)
{
  return OpMul<SparseNonZeroPattern<TM>, SparseVector<TV>>( A, x );
}

} //namespace SLA
} //namespace SANS


#endif //DENSEMATRIX_MUL_H
