// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXD_ADD_H
#define MATRIXD_ADD_H

#include "tools/SANSnumerics.h"
#include "MatrixD_Type.h"

namespace numpack 
{
namespace DLA
{

// Lazy expressions based on recursive function calls

// Addition and Subtraction
// No temporary variables are needed for these operations

template<class ExprL, class ExprR, bool useRF>
class OpAddD : public MatrixDType< OpAddD<ExprL,ExprR, useRF>, useRF >
{
public:
  typedef typename ExprL::node_type node_type;

  OpAddD(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR)
  {
    SANS_ASSERT( eL.m() == eR.m() );
    SANS_ASSERT( eL.n() == eR.n() );
  }

  template< class T >
  inline void value(const Real sgn, MatrixDView<T>& res) const
  {
    eL.value(sgn, res);
    eR.plus(sgn, res);
  }
  template< class T >
  inline void plus(const Real sgn, MatrixDView<T>& res) const
  {
    eL.plus(sgn, res);
    eR.plus(sgn, res);
  }

  template< class MatrixL >
  inline void value(const Real sgn, MatrixDTuple<MatrixL>& res) const
  {
    eL.value(sgn, res);
    eR.plus(sgn, res);
  }
  template< class MatrixL >
  inline void plus(const Real sgn, MatrixDTuple<MatrixL>& res) const
  {
    eL.plus(sgn, res);
    eR.plus(sgn, res);
  }

  //Element-wise expression
  inline node_type value(const int& i) const
  {
    return eL.value(i) + eR.value(i);
  }

  inline const OpAddD&
  operator+() const { return *this; }
  int m() const { return eL.m(); }
  int n() const { return eL.n(); }
  int size() const { return m()*n(); }
private:
  const ExprL& eL;
  const ExprR& eR;
};


template<class ExprL, class ExprR, bool useRF>
class OpSubD : public MatrixDType< OpSubD<ExprL, ExprR, useRF>, useRF >
{
public:
  typedef typename ExprL::node_type node_type;

  OpSubD(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR)
  {
    SANS_ASSERT( eL.m() == eR.m() );
    SANS_ASSERT( eL.n() == eR.n() );
  }

  template< class T >
  inline void value(const Real sgn, MatrixDView<T>& res) const
  {
    eL.value(sgn, res);
    eR.plus(-sgn, res);
  }
  template< class T >
  inline void plus(const Real sgn, MatrixDView<T>& res) const
  {
    eL.plus(sgn, res);
    eR.plus(-sgn, res);
  }

  template< class MatrixL >
  inline void value(const Real sgn, MatrixDTuple<MatrixL>& res) const
  {
    eL.value(sgn, res);
    eR.plus(-sgn, res);
  }
  template< class MatrixL >
  inline void plus(const Real sgn, MatrixDTuple<MatrixL>& res) const
  {
    eL.plus(sgn, res);
    eR.plus(-sgn, res);
  }

  //Element-wise expression
  inline node_type value(const int& i) const
  {
    return eL.value(i) - eR.value(i);
  }

  inline const OpSubD&
  operator+() const { return *this; }
  int m() const { return eL.m(); }
  int n() const { return eL.n(); }
  int size() const { return m()*n(); }
private:
  const ExprL& eL;
  const ExprR& eR;
};

//Generator functions to generate an addition/subtraction operation
template<class ExprL, bool useRFL, class ExprR, bool useRFR>
inline OpAddD<ExprL, ExprR, useRFL || useRFR>
operator+(const MatrixDType<ExprL, useRFL>& eL, const MatrixDType<ExprR, useRFR>& eR)
{
  return OpAddD<ExprL, ExprR, useRFL || useRFR>( eL.cast(), eR.cast() );
}


template<class ExprL, bool useRFL, class ExprR, bool useRFR>
inline OpSubD<ExprL, ExprR, useRFL || useRFR>
operator-(const MatrixDType<ExprL, useRFL>& eL, const MatrixDType<ExprR, useRFR>& eR)
{
  return OpSubD<ExprL, ExprR, useRFL || useRFR>( eL.cast(), eR.cast() );
}


} //namespace DLA
} //namespace numpack 


/* This is dissabled until we decide that it's a good idea to add/subtract scalar quantities with a matrix
namespace numpack 
{
namespace DLA
{

//Addition and Subtraction with scalar quantities

template<class Expr>
class OpScalar : public MatrixDType< OpScalar<Expr> >
{
public:
  typedef typename Expr::node_type node_type;

  OpScalar(const Expr& e, const Real esgn, const Real s) : e(e), esgn(esgn), s(s) {}

  template< class T >
  inline void value(const Real sgn, MatrixDView<T>& res) const
  {
    e.value(esgn*sgn, res);
    res.value() += sgn*s;
  }
  template< class T >
  inline void plus(const Real sgn, MatrixDView<T>& res) const
  {
    e.plus(esgn*sgn, res);
    res.value() += sgn*s;
  }

  inline const OpScalar&
  operator+() const { return *this; }
  int m() const { return e.m(); }
  int n() const { return e.n(); }
  int size() const { return m()*n(); }
private:
  const Expr& e;
  const Real esgn;
  const Real s;
};

}
}


template<class Expr>
inline SurrealExp::OpScalar<Expr>
operator+(const MatrixDType<Expr>& e, const Real& s)
{
  return SurrealExp::OpScalar<Expr>( e.cast(), 1, s );
}

template<class Expr>
inline SurrealExp::OpScalar<Expr>
operator+(const Real& s, const MatrixDType<Expr>& e)
{
  return SurrealExp::OpScalar<Expr>( e.cast(), 1, s );
}

template<class Expr>
inline SurrealExp::OpScalar<Expr>
operator-(const MatrixDType<Expr>& e, const Real& s)
{
  return SurrealExp::OpScalar<Expr>( e.cast(), 1, -s );
}

template<class Expr>
inline SurrealExp::OpScalar<Expr>
operator-(const Real& s, const MatrixDType<Expr>& e)
{
  return SurrealExp::OpScalar<Expr>( e.cast(), -1, s );
}
*/

#endif //MATRIXD_ADD_H
