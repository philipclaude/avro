// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef sparse_ADD_H
#define sparse_ADD_H

#include "tools/SANSnumerics.h"
#include "sparse_Type.h"

#include "SparseVector.h"

namespace numpack 
{
namespace SLA
{

// Lazy expressions based on recursive function calls

// Addition and Subtraction
// No temporary variables are needed for these operations
template<class ExprL, class ExprR, bool useRF>
class OpAdd;

template<class ExprL, class ExprR>
class OpAdd<ExprL, ExprR, true> : public sparseType< OpAdd<ExprL, ExprR, true>, true >
{
public:
  typedef typename ExprL::Ttype Ttype;

  OpAdd(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR) { SANS_ASSERT( eL.m() == eR.m() ); }

  template< class T >
  inline void value(const Real sgn, SparseVector<T>& res) const
  {
    eL.value(sgn, res);
    eR.plus(sgn, res);
  }
  template< class T >
  inline void plus(const Real sgn, SparseVector<T>& res) const
  {
    eL.plus(sgn, res);
    eR.plus(sgn, res);
  }
/*
  template< class MatrixL >
  inline void value(const Real sgn, sparseTuple<MatrixL>& res) const
  {
    Ll.value(sgn, res);
    Rr.plus(sgn, res);
  }
  template< class MatrixL >
  inline void plus(const Real sgn, sparseTuple<MatrixL>& res) const
  {
    Ll.plus(sgn, res);
    Rr.plus(sgn, res);
  }
*/

  inline const OpAdd&
  operator+() const { return *this; }
  int m() const { return eL.m(); }
private:
  const ExprL& eL;
  const ExprR& eR;
};

template<class ExprL, class ExprR>
class OpAdd<ExprL, ExprR, false> : public sparseType< OpAdd<ExprL, ExprR, false>, false >
{
public:
  typedef typename ExprL::Ttype Ttype;

  OpAdd(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR) { SANS_ASSERT( eL.m() == eR.m() ); }

  //Element-wise expression
  inline Ttype operator[](const int& i) const
  {
    return eL[i] + eR[i];
  }

  inline const OpAdd&
  operator+() const { return *this; }
  int m() const { return eL.m(); }
private:
  const ExprL& eL;
  const ExprR& eR;
};

template<class ExprL, class ExprR, bool useRF>
class OpSub;

template<class ExprL, class ExprR>
class OpSub<ExprL, ExprR, true> : public sparseType< OpSub<ExprL, ExprR, true>, true >
{
public:
  typedef typename ExprL::Ttype Ttype;

  OpSub(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR) { SANS_ASSERT( eL.m() == eR.m() ); }

  template< class T >
  inline void value(const Real sgn, SparseVector<T>& res) const
  {
    eL.value(sgn, res);
    eR.plus(-sgn, res);
  }
  template< class T >
  inline void plus(const Real sgn, SparseVector<T>& res) const
  {
    eL.plus(sgn, res);
    eR.plus(-sgn, res);
  }
/*
  template< class MatrixL >
  inline void value(const Real sgn, sparseTuple<MatrixL>& res) const
  {
    Ll.value(sgn, res);
    Rr.plus(-sgn, res);
  }
  template< class MatrixL >
  inline void plus(const Real sgn, sparseTuple<MatrixL>& res) const
  {
    Ll.plus(sgn, res);
    Rr.plus(-sgn, res);
  }
*/

  inline const OpSub&
  operator+() const { return *this; }
  int m() const { return eL.m(); }
private:
  const ExprL& eL;
  const ExprR& eR;
};


template<class ExprL, class ExprR>
class OpSub<ExprL, ExprR, false> : public sparseType< OpSub<ExprL, ExprR, false>, false >
{
public:
  typedef typename ExprL::Ttype Ttype;

  OpSub(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR) { SANS_ASSERT( eL.m() == eR.m() ); }

  //Element-wise expression
  inline Ttype operator[](const int& i) const
  {
    return eL[i] - eR[i];
  }

  inline const OpSub&
  operator+() const { return *this; }
  int m() const { return eL.m(); }
private:
  const ExprL& eL;
  const ExprR& eR;
};

//Generator functions to generate an addition/subtraction operation
template<class ExprL, bool useRFL, class ExprR, bool useRFR>
inline OpAdd<ExprL, ExprR, useRFL || useRFR>
operator+(const sparseType<ExprL, useRFL>& eL, const sparseType<ExprR, useRFR>& eR)
{
  return OpAdd<ExprL, ExprR, useRFL || useRFR>( eL.cast(), eR.cast() );
}

template<class ExprL, bool useRFL, class ExprR, bool useRFR>
inline OpSub<ExprL, ExprR, useRFL || useRFR>
operator-(const sparseType<ExprL, useRFL>& eL, const sparseType<ExprR, useRFR>& eR)
{
  return OpSub<ExprL, ExprR, useRFL || useRFR>( eL.cast(), eR.cast() );
}


} //namespace SLA
} //namespace numpack 

#endif //sparse_ADD_H
