// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef BLOCKLINALG_ADD_H
#define BLOCKLINALG_ADD_H

#include "tools/SANSnumerics.h"
#include "tools/SANSException.h"

#include "BlockLinAlg_Type.h"

namespace numpack 
{
namespace BLA
{
// Addition and Subtraction of BlockLinAlgType (vector, matrix or expression)
//
// Lazy expressions based on recursive function calls
// No temporary variables are needed for these operations
//
//---------------------------------------------------------------------------//
template<class ExprL, class ExprR>
class OpAdd : public BlockLinAlgType< OpAdd<ExprL, ExprR> >
{
public:
  OpAdd(const ExprL& eL, const ExprR& eR) : eL_(eL), eR_(eR)
  {
    SANS_ASSERT(eL_.m() == eR_.m());
    SANS_ASSERT(eL_.n() == eR_.n());
  }

  template< class Vector >
  inline void value(const Real sgn, BlockVectorType<Vector>& res) const
  {
    eL_.value(sgn, res.cast());
    eR_.plus(sgn, res.cast());
  }

  template< class Vector >
  inline void plus(const Real sgn, BlockVectorType<Vector>& res) const
  {
    eL_.plus(sgn, res.cast());
    eR_.plus(sgn, res.cast());
  }

  template< class Matrix >
  inline void value(const Real sgn, BlockMatrixType<Matrix>& res) const
  {
    eL_.value(sgn, res.cast());
    eR_.plus(sgn, res.cast());
  }

  template< class Matrix >
  inline void plus(const Real sgn, BlockMatrixType<Matrix>& res) const
  {
    eL_.plus(sgn, res.cast());
    eR_.plus(sgn, res.cast());
  }

  inline const OpAdd& operator+() const { return *this; }

  int m() const { return eL_.m(); }
  int n() const { return eL_.n(); }

private:
  const ExprL& eL_;
  const ExprR& eR_;
};

//---------------------------------------------------------------------------//
template<class ExprL, class ExprR>
class OpSub : public BlockLinAlgType< OpSub<ExprL, ExprR> >
{
public:
  OpSub(const ExprL& eL, const ExprR& eR) : eL_(eL), eR_(eR)
  {
    SANS_ASSERT(eL_.m() == eR_.m());
    SANS_ASSERT(eL_.n() == eR_.n());
  }

  template< class Vector >
  inline void value(const Real sgn, BlockVectorType<Vector>& res) const
  {
    eL_.value(sgn, res.cast());
    eR_.plus(-sgn, res.cast());
  }
  template< class Vector >
  inline void plus(const Real sgn, BlockVectorType<Vector>& res) const
  {
    eL_.plus(sgn, res.cast());
    eR_.plus(-sgn, res.cast());
  }

  template< class Matrix >
  inline void value(const Real sgn, BlockMatrixType<Matrix>& res) const
  {
    eL_.value(sgn, res.cast());
    eR_.plus(-sgn, res.cast());
  }

  template< class Matrix >
  inline void plus(const Real sgn, BlockMatrixType<Matrix>& res) const
  {
    eL_.plus(sgn, res.cast());
    eR_.plus(-sgn, res.cast());
  }

  inline const OpSub& operator+() const { return *this; }

  int m() const { return eL_.m(); }
  int n() const { return eL_.n(); }

private:
  const ExprL& eL_;
  const ExprR& eR_;
};

//---------------------------------------------------------------------------//
// Overloaded operators for generating a representation of addition/subtraction operation
template<class ExprL, class ExprR>
inline OpAdd<ExprL, ExprR>
operator+(const BlockLinAlgType<ExprL>& eL, const BlockLinAlgType<ExprR>& eR)
{
  return OpAdd<ExprL, ExprR>(eL.cast(), eR.cast());
}

template<class ExprL, class ExprR>
inline OpSub<ExprL, ExprR>
operator-(const BlockLinAlgType<ExprL>& eL, const BlockLinAlgType<ExprR>& eR)
{
  return OpSub<ExprL, ExprR>(eL.cast(), eR.cast());
}

} //namespace BLA
} //namespace numpack 

#endif //BLOCKLINALG_ADD_H
