// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXD_MUL_H
#define MATRIXD_MUL_H

// Surreals are registered with std::is_arithmetic.
// not sure how to register them with std::is_arithmetic
//#include <boost/type_traits/is_arithmetic.hpp>

#include <type_traits>

#include "MatrixD_Type.h"
#include "MatMul/MatMul_BLAS.h"
#include "MatMul/MatMul_Native.h"


namespace numpack
{
namespace SLA
{
template<class TM> class SparseMatrix_CRS;
template<class TV> class SparseVector;
}

namespace DLA
{

//Multiplication with MatrixD

//-----------------------------------------------------------------------------
// OpMulD_impl implements the actual multiplication algorithm
template<class TL, class TR, class T>
class OpMulD_impl
{
public:
  static inline void value(const MatrixDView<TL>& ML, const MatrixDView<TR>& MR, const Real sgn, MatrixDView<T>& res )
  {
    MatMul_Native<TL,TR,T>::value(ML, MR, sgn, res);
  }

  static inline void plus(const MatrixDView<TL>& ML, const MatrixDView<TR>& MR, const Real sgn, MatrixDView<T>& res )
  {
    MatMul_Native<TL,TR,T>::plus(ML, MR, sgn, res);
  }

};

template<class TM, class TV, class T>
class OpMulD_impl< SLA::SparseMatrix_CRS<TM>, SLA::SparseVector<TV>, SLA::SparseVector<T> >
{
public:
  static inline void value(const MatrixDView< SLA::SparseMatrix_CRS<TM> >& A,
                           const MatrixDView< SLA::SparseVector<TV> >& x,
                           const Real sgn,
                           MatrixDView< SLA::SparseVector<T> >& res )
  {
    res = 0;
    plus(A, x, sgn, res);
  }

  static inline void plus(const MatrixDView< SLA::SparseMatrix_CRS<TM> >& A,
                          const MatrixDView< SLA::SparseVector<TV> >& x,
                          const Real sgn,
                          MatrixDView< SLA::SparseVector<T> >& res )
  {
    //No need to do any optimization here as we are likely working with large sparse matrices anyways
    SANS_ASSERT( A.m() == res.m() );
    SANS_ASSERT( A.n() == x.m() );
    SANS_ASSERT( x.n() == 1 && res.n() == 1 );

    for (int i = 0; i < A.m(); i++)
      for (int j = 0; j < A.n(); j++)
        res(i,0) += sgn*(A(i,j)*x(j,0));
  }
};

#ifdef DLA_BLAS
//-----------------------------------------------------------------------------
template<>
class OpMulD_impl<Real, Real, Real>
{
public:
  static void value(const MatrixDView<Real>& ML, const MatrixDView<Real>& MR, const Real sgn, MatrixDView<Real>& res );

  static void plus(const MatrixDView<Real>& ML, const MatrixDView<Real>& MR, const Real sgn, MatrixDView<Real>& res );
};
#endif

//-----------------------------------------------------------------------------
// This is a general multiplication of two expressions, i.e.
// M = (A + B)*(C + D)
// The expressions must be evaluated into temporary variables in order to evaluate
// the multiplication
template<class ExprL, class ExprR>
class OpMulD : public MatrixDType< OpMulD<ExprL, ExprR>, true >
{
public:
  typedef typename ExprL::node_type TL;
  typedef typename ExprR::node_type TR;
  typedef TR node_type;

  OpMulD(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR)
  {
    SANS_ASSERT( eL.n() == eR.m() );
  }

  template<class T>
  inline void value(const Real sgn, MatrixDView<T>& res) const
  {
    MatrixD<TL> ML( eL );
    MatrixD<TR> MR( eR );

    OpMulD_impl<TL, TR, T>::value( ML, MR, sgn, res );
  }
  template<class T>
  inline void plus(const Real sgn, MatrixDView<T>& res) const
  {
    MatrixD<TL> ML( eL );
    MatrixD<TR> MR( eR );

    OpMulD_impl<TL, TR, T>::plus( ML, MR, sgn, res );
  }

  template< class MatrixL >
  inline void value(const Real sgn, MatrixDTuple<MatrixL>& tuple) const
  {
    MatrixD<TL> ML( eL );
    MatrixD<TR> MR( eR );
    MatrixD<TR> res( tuple.m(), tuple.n() );

    OpMulD_impl<TL, TR, TR>::value( ML, MR, sgn, res );
    tuple = res;
  }
  template< class MatrixL >
  inline void plus(const Real sgn, MatrixDTuple<MatrixL>& tuple) const
  {
    MatrixD<TL> ML( eL );
    MatrixD<TR> MR( eR );
    MatrixD<TR> res( tuple );

    OpMulD_impl<TL, TR, TR>::plus( ML, MR, sgn, res );
    tuple = res;
  }

  operator TL() const
  {
    SANS_ASSERT( m() == 1 );
    SANS_ASSERT( n() == 1 );
    MatrixD<TL> M(*this);
    return M(0,0);
  }

  inline const OpMulD&
  operator+() const { return *this; }
  int m() const { return eL.m(); }
  int n() const { return eR.n(); }
  int size() const { return m()*n(); }
private:
  const ExprL& eL;
  const ExprR& eR;
};

//-----------------------------------------------------------------------------
// This specialization is for multiplication between a matrix and an expression, i.e.
// M = A*(B + C)
// only the expression needs to be evaluated in a temporary variable
template<class TL, class ExprR>
class OpMulD< MatrixDView<TL>, ExprR > : public MatrixDType< OpMulD< MatrixDView<TL>, ExprR >, true >
{
public:
  typedef typename ExprR::node_type TR;
  typedef TR node_type;

  OpMulD(const MatrixDView<TL>& ML, const ExprR& eR) :  ML(ML), eR(eR)
  {
    SANS_ASSERT( ML.n() == eR.m() );
  }

  template<class T>
  inline void value(const Real sgn, MatrixDView<T>& res) const
  {
    MatrixD<TR> MR( eR );
    OpMulD_impl<TL,TR,T>::value( ML, MR, sgn, res );
  }
  template<class T>
  inline void plus(const Real sgn, MatrixDView<T>& res) const
  {
    MatrixD<TR> MR( eR );
    OpMulD_impl<TL,TR,T>::plus( ML, MR, sgn, res );
  }

  template< class MatrixR >
  inline void value(const Real sgn, MatrixDTuple<MatrixR>& tuple) const
  {
    MatrixD<TR> MR( eR );
    MatrixD<TR> res( tuple.m(), tuple.n() );

    OpMulD_impl<TL,TR,TR>::value( ML, MR, sgn, res );
    tuple = res;
  }
  template< class MatrixR >
  inline void plus(const Real sgn, MatrixDTuple<MatrixR>& tuple) const
  {
    MatrixD<TR> MR( eR );
    MatrixD<TR> res( tuple );

    OpMulD_impl<TL,TR,TR>::plus( ML, MR, sgn, res );
    tuple = res;
  }

  operator TL() const
  {
    SANS_ASSERT( m() == 1 );
    SANS_ASSERT( n() == 1 );
    MatrixD<TL> M(*this);
    return M(0,0);
  }

  inline const OpMulD&
  operator+() const { return *this; }
  int m() const { return ML.m(); }
  int n() const { return eR.n(); }
  int size() const { return m()*n(); }
private:
  const MatrixDView<TL>& ML;
  const ExprR& eR;
};

//-----------------------------------------------------------------------------
// This specialization is for multiplication between an expression and a matrix, i.e.
// M = (A + B)*C
// only the expression needs to be evaluated in a temporary variable
template<class ExprL, class TR>
class OpMulD< ExprL, MatrixDView<TR> > : public MatrixDType< OpMulD< ExprL, MatrixDView<TR> >, true >
{
public:
  typedef typename ExprL::node_type TL;
  typedef TR node_type;

  OpMulD(const ExprL& eL, const MatrixDView<TR>& MR) : eL(eL), MR(MR)
  {
    SANS_ASSERT( eL.n() == MR.m() );
  }

  template<class T>
  inline void value(const Real sgn, MatrixDView<T>& res) const
  {
    MatrixD<TL> ML( eL );
    OpMulD_impl<TL, TR, T>::value( ML, MR, sgn, res );
  }
  template<class T>
  inline void plus(const Real sgn, MatrixDView<T>& res) const
  {
    MatrixD<TL> ML( eL );
    OpMulD_impl<TL, TR, T>::plus( ML, MR, sgn, res );
  }

  /*
  //removed by pcaplan sept 21, 2019 (doesn't seem to be used in matrixs_mul.h)
  template< class MatrixR >
  inline void value(const Real sgn, MatrixDTuple<MatrixR>& tuple) const
  {
    MatrixD<TL> ML( eL );
    MatrixD<TR> res( tuple.m(), tuple.n() );

    OpMulD_impl<TL, TR, TR>::value( ML, MR, sgn, res );
    tuple = res;
  }
  template< class MatrixR >
  inline void plus(const Real sgn, MatrixDTuple<MatrixR>& tuple) const
  {
    MatrixD<TL> ML( eL );
    MatrixD<TR> res( tuple );

    OpMulD_impl<TL, TR, TR>::plus( ML, MR, sgn, res );
    tuple = res;
  }
  */

  operator TL() const
  {
    SANS_ASSERT( m() == 1 );
    SANS_ASSERT( n() == 1 );
    MatrixD<TL> M(*this);
    return M(0,0);
  }

  inline const OpMulD&
  operator+() const { return *this; }
  int m() const { return eL.m(); }
  int n() const { return MR.n(); }
  int size() const { return m()*n(); }
private:
  const ExprL& eL;
  const MatrixDView<TR>& MR;
};

//-----------------------------------------------------------------------------
// This specialization is for multiplication between two matrices, i.e.
// M = A*B
// no temporary variables are required
template<class TL, class TR>
class OpMulD< MatrixDView<TL>, MatrixDView<TR> > : public MatrixDType< OpMulD< MatrixDView<TL>, MatrixDView<TR> >, true >
{
public:
  typedef TR node_type;

  OpMulD(const MatrixDView<TL>& ML, const MatrixDView<TR>& MR) : ML(ML), MR(MR)
  {
    SANS_ASSERT( ML.n() == MR.m() );
  }

  template<class T>
  inline void value(const Real sgn, MatrixDView<T>& res) const
  {
    OpMulD_impl<TL, TR, T>::value( ML, MR, sgn, res );
  }
  template<class T>
  inline void plus(const Real sgn, MatrixDView<T>& res) const
  {
    OpMulD_impl<TL, TR, T>::plus( ML, MR, sgn, res );
  }

  //A temporary variable cannot be avoided when res is a MatrixTuple
  template< class MatrixR >
  inline void value(const Real sgn, MatrixDTuple<MatrixR >& tuple) const
  {
    MatrixD<TR> res( tuple.m(), tuple.n() );
    OpMulD_impl<TL, TR, TR>::value( ML, MR, sgn, res );
    tuple = res;
  }
  template< class MatrixR >
  inline void plus(const Real sgn, MatrixDTuple<MatrixR>& tuple) const
  {
    MatrixD<TR> res( tuple );
    OpMulD_impl<TL, TR, TR>::plus( ML, MR, sgn, res );
    tuple = res;
  }

  operator TL() const
  {
    SANS_ASSERT( m() == 1 );
    SANS_ASSERT( n() == 1 );
    MatrixD<TL> M(*this);
    return M(0,0);
  }

  inline const OpMulD&
  operator+() const { return *this; }
  int m() const { return ML.m(); }
  int n() const { return MR.n(); }
  int size() const { return m()*n(); }
private:
  const MatrixDView<TL>& ML;
  const MatrixDView<TR>& MR;
};


//-----------------------------------------------------------------------------
// Operator for generating an OpMulD representation of matrix multiplication including expressions
template<class ExprL, bool useRFL, class ExprR, bool useRFR>
inline OpMulD<ExprL, ExprR>
operator*(const MatrixDType<ExprL, useRFL>& eL, const MatrixDType<ExprR, useRFR>& eR)
{
  return OpMulD<ExprL, ExprR>( eL.cast(), eR.cast() );
}


// Multiplication between a MatrixD expression and a scalar quantity


//-----------------------------------------------------------------------------
// This represents multiplication between a matrix or matrix expression and a scalar, i.e.
// M = s*(A + B)
// or
// M = (A + B)*s
// where s is a scalar quantity. The class also represents division with s by storing 1/s.
template<class Expr, bool useRF>
class OpMulDScalar : public MatrixDType< OpMulDScalar<Expr, useRF>, useRF >
{
public:
  typedef typename Expr::node_type node_type;

  OpMulDScalar(const Expr& e, const Real s) : e(e), s(s) {}

  template<class T>
  inline void value(const Real sgn, MatrixDView<T>& res) const
  {
    e.value(s*sgn, res);
  }
  template<class T>
  inline void plus(const Real sgn, MatrixDView<T>& res) const
  {
    e.plus(s*sgn, res);
  }

  //Element-wise expression
  inline node_type value(const int& i) const
  {
    return s*e.value(i);
  }

  inline const OpMulDScalar&
  operator+() const { return *this; }
  int m() const { return e.m(); }
  int n() const { return e.n(); }
  int size() const { return m()*n(); }

  const Expr& e;
  const Real s;
};


//=============================================================================
// Overloaded operators to represent multiplication between a scalar and a matrix expression
template<class Expr, bool useRF, typename T>
inline typename std::enable_if< std::is_arithmetic<T>::value, OpMulDScalar<Expr, useRF> >::type
operator*(const MatrixDType<Expr, useRF>& e, const T& s)
{
  return OpMulDScalar<Expr, useRF>( e.cast(), s );
}

template<class Expr, bool useRF, typename T>
inline typename std::enable_if< std::is_arithmetic<T>::value, OpMulDScalar<Expr, useRF> >::type
operator/(const MatrixDType<Expr, useRF>& e, const T& s)
{
  return OpMulDScalar<Expr, useRF>( e.cast(), Real(1)/s );
}

template<class Expr, bool useRF, typename T>
inline typename std::enable_if< std::is_arithmetic<T>::value, OpMulDScalar<Expr, useRF> >::type
operator*(const T& s, const MatrixDType<Expr, useRF>& e)
{
  return OpMulDScalar<Expr, useRF>( e.cast(), s );
}

//=============================================================================
// This is a special case when multiplies scalars are multiplying from two sides, i.e. B = 2*A*2;
// This reduces the complexity of the expression tree and hence reduces code bloat
template<class Expr, bool useRF, typename T>
inline typename std::enable_if< std::is_arithmetic<T>::value, OpMulDScalar<Expr, useRF> >::type
operator*(const OpMulDScalar<Expr, useRF>& MulScal, const T& s)
{
  return OpMulDScalar<Expr, useRF>( MulScal.e, MulScal.s*s );
}

template<class Expr, bool useRF, typename T>
inline typename std::enable_if< std::is_arithmetic<T>::value, OpMulDScalar<Expr, useRF> >::type
operator/(const OpMulDScalar<Expr, useRF>& MulScal, const T& s)
{
  return OpMulDScalar<Expr, useRF>( MulScal.e, MulScal.s/s );
}

template<class Expr, bool useRF, typename T>
inline typename std::enable_if< std::is_arithmetic<T>::value, OpMulDScalar<Expr, useRF> >::type
operator*(const T& s, const OpMulDScalar<Expr, useRF>& MulScal)
{
  return OpMulDScalar<Expr, useRF>( MulScal.e, MulScal.s*s );
}


template<class ExprL, class ExprR>
class OpMulDFactor;

//-----------------------------------------------------------------------------
// This specialization is for multiplication that involve a scalar and translates
// (s*A)*x into s*(A*x)
// no temporary variables are required
template<class ExprL, bool useRFL, class ExprR>
class OpMulDFactor< OpMulDScalar<ExprL, useRFL>, ExprR >
  : public MatrixDType< OpMulDFactor< OpMulDScalar<ExprL, useRFL>, ExprR >, true >
{
public:
  typedef typename ExprL::node_type node_type;

  OpMulDFactor(const OpMulDScalar<ExprL, useRFL>& eL, const ExprR& eR) : eL(eL), eR(eR) {}

  template<class TR>
  inline void value(const Real sgn, MatrixDView<TR>& res) const
  {
    res = eL.s*sgn*(eL.e*eR);
  }
  template<class TR>
  inline void plus(const Real sgn, MatrixDView<TR>& res) const
  {
    res += eL.s*sgn*(eL.e*eR);
  }

  inline const OpMulDFactor&
  operator+() const { return *this; }
  int m() const { return eL.m(); }
  int n() const { return eR.n(); }
  int size() const { return m()*n(); }
private:
  const OpMulDScalar<ExprL, useRFL>& eL;
  const ExprR& eR;
};

//-----------------------------------------------------------------------------
// This specialization is for multiplication that involve a scalar and translates
// A*(x*s) into (A*x)*s
// no temporary variables are required
template<class ExprL, class ExprR, bool useRFR>
class OpMulDFactor< ExprL, OpMulDScalar<ExprR, useRFR> >
  : public MatrixDType< OpMulDFactor< ExprL, OpMulDScalar<ExprR, useRFR> >, true >
{
public:
  typedef typename ExprL::node_type node_type;

  OpMulDFactor(const ExprL& eL, const OpMulDScalar<ExprR, useRFR>& eR) : eL(eL), eR(eR) {}

  template<class TV>
  inline void value(const Real sgn, MatrixDView<TV>& res) const
  {
    res = eR.s*sgn*(eL*eR.e);
  }
  template<class TV>
  inline void plus(const Real sgn, MatrixDView<TV>& res) const
  {
    res += eR.s*sgn*(eL*eR.e);
  }

  inline const OpMulDFactor&
  operator+() const { return *this; }
  int m() const { return eL.m(); }
  int n() const { return eR.n(); }
  int size() const { return m()*n(); }
private:
  const ExprL& eL;
  const OpMulDScalar<ExprR, useRFR>& eR;
};

//-----------------------------------------------------------------------------
// Operator for generating an OpMulD representation of matrix multiplication including expressions
template<class ExprL, bool useRFL, class ExprR, bool useRFR>
inline OpMulDFactor< OpMulDScalar<ExprL, useRFL>, ExprR>
operator*(const OpMulDScalar<ExprL, useRFL>& eL, const MatrixDType<ExprR, useRFR>& eR)
{
  return OpMulDFactor< OpMulDScalar<ExprL, useRFL>, ExprR>( eL, eR.cast() );
}

template<class ExprL, bool useRFL, class ExprR, bool useRFR>
inline OpMulDFactor<ExprL, OpMulDScalar<ExprR, useRFR> >
operator*(const MatrixDType<ExprL, useRFL>& eL, const OpMulDScalar<ExprR, useRFR>& eR)
{
  return OpMulDFactor<ExprL, OpMulDScalar<ExprR, useRFR> >( eL.cast(), eR );
}

} //namespace DLA
} //namespace numpack


#endif //MATRIXD_MUL_H
