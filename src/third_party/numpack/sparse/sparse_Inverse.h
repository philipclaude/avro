// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef sparse_INVERSE_H
#define sparse_INVERSE_H

#include <memory>

#include "numpack/dense/dynamic/MatrixD_Type.h"
#include "numpack/block/block_Type.h"
#include "sparse_Type.h"

#include "numpack/VectorType.h"

#include "LinearSolverBase.h"
#include "SparseVector.h"
#include "sparse_Mul.h"

namespace numpack 
{
namespace SLA
{
#if 0
//=============================================================================
//Represents the multiplication of a matrix inverse with a vector, i.e. Solver.Inverse(A)*b;
template< class Matrix_type >
class SparseBackSolverMul : public sparseType< SparseBackSolverMul< Matrix_type >, true >,
                            public DLA::MatrixDType< SparseBackSolverMul< Matrix_type >, true >,
                            public BLA::blockType< SparseBackSolverMul< Matrix_type > >
{
public:

  typedef typename VectorType<Matrix_type>::Viewtype SparseVectorView_type;

  SparseBackSolverMul( const SparseBackSolver<Matrix_type>& InvA, const SparseVectorView_type& b, const Real s = 1 )
    : InvA_(InvA), b_(b), s_(s) {}

  //This allows for the general expression x = s*Inverse( A )*b
  inline void value(const Real sgn, SparseVectorView_type& x) const
  {
    InvA_.Solver.backsolve(b_, s_*sgn, x);
  }

  inline const SparseBackSolverMul&
  operator+() const { return *this; }
  int m() const { return InvA_.n(); }
  int n() const { return b_.n(); }
protected:
  const SparseBackSolver<Matrix_type>& InvA_; //The sparse matrix solver
  const SparseVectorView_type& b_;         //The sparse vector multiplying the matrix inverse
  Real s_;                                 //A scalar quantity multiplying the inverse operation
};

//=============================================================================
//Operator for multiplication between an inverse-matrix and a matrix, i.e. Solver.Inverse(A)*b;
template< class Matrix_type >
inline SparseBackSolverMul< Matrix_type >
operator*(const SparseBackSolver< Matrix_type >& InvA,
          const typename VectorType<Matrix_type>::type& b )
{
  return SparseBackSolverMul< Matrix_type >( InvA, b );
}

//=============================================================================
//Operator for multiplication between a scalar, an inverse-matrix, and a matrix, i.e. 2*Solver.Inverse(A)*b;
template< class Matrix_type, bool useRF >
inline SparseBackSolverMul< Matrix_type >
operator*(const OpMulScalar< SparseBackSolver< Matrix_type >, useRF >& ScalMulInvA,
          const typename VectorType<Matrix_type>::type& b )
{
  return SparseBackSolverMul< Matrix_type >( ScalMulInvA.e, b, ScalMulInvA.s);
}
#endif
} //namespace SLA
} //namespace numpack 

#endif //sparse_INVERSE_H
