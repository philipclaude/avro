// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef ILU0_H
#define ILU0_H

#include "Python/PyDict.h"

#include "tools/noncopyable.h"

#include "LinearAlgebra/AlgebraicEquationSetBase.h"
#include "LinearAlgebra/SparseLinAlg/LinearSolverBase.h"
#include "LinearAlgebra/SparseLinAlg/SparseLinAlg_Inverse.h"
#include "LinearAlgebra/SparseLinAlg/SparseMatrix_Diag.h"

#include "tools/minmax.h"

#include <vector>

namespace SANS
{
namespace SLA
{


//=============================================================================
struct ILU0Param : noncopyable
{
  static void checkInputs(PyDict d);
  static ILU0Param params;
};


//=============================================================================
//
// Incomplete Lower-Upper Factorization with 0 fill in
//
// The matrix is decomposed into lower, L, upper, U, and remainder, R, matrices as
//
//   A = LU + R
//
// The ILU preconditioner approximates A as
//
// M = LU
//
template< class Matrix_type >
class ILU0 : public LinearSolverBase< Matrix_type >
{
public:
  typedef LinearSolverBase< Matrix_type > Base_type;

  typedef typename Base_type::SparseVector_type SparseVector_type;
  typedef typename Base_type::SparseVectorView_type SparseVectorView_type;
  typedef typename NonZeroPatternType<Matrix_type>::type SystemNonZeroPattern;

  typedef typename Matrix_type::Ttype TM;
  typedef typename SparseVector_type::Ttype TV;

//-----------------------------------------------------------------------------
  ILU0( const PyDict& d, AlgebraicEquationSetBase<Matrix_type>& f, LinearSystemSolve solve = RegularSolve) :
    Base_type(solve), f_(f)
  {
    init();
  }
  virtual ~ILU0() { delete A_; }

//-----------------------------------------------------------------------------
  virtual void factorize() override;

//-----------------------------------------------------------------------------
  virtual LinearSolveStatus backsolve( const SparseVectorView_type& b, SparseVectorView_type& x ) const override;

protected:
  void init();

  using Base_type::A_;
  using Base_type::transpose_;
  AlgebraicEquationSetBase<Matrix_type>& f_;
  SparseMatrix_CRS<TM> LU_;
};



//=============================================================================
template< class Matrix_type >
class ILU0< DLA::MatrixD<Matrix_type> > : public LinearSolverBase< DLA::MatrixD<Matrix_type> >
{
public:
  typedef LinearSolverBase< DLA::MatrixD<Matrix_type> > Base_type;

  typedef typename Base_type::SparseVector_type SparseVector_type;
  typedef typename Base_type::SparseVectorView_type SparseVectorView_type;
  typedef typename NonZeroPatternType<DLA::MatrixD<Matrix_type>>::type SystemNonZeroPattern;

  typedef typename Matrix_type::Ttype TM;
  typedef typename SparseVector_type::node_type::Ttype TV;

//-----------------------------------------------------------------------------
  ILU0( const PyDict& d, AlgebraicEquationSetBase<DLA::MatrixD<Matrix_type>>& f, LinearSystemSolve solve = RegularSolve ) :
    Base_type(solve), f_(f)
  {
    init();
  }
  virtual ~ILU0() {}

//-----------------------------------------------------------------------------
  virtual void factorize() override;

//-----------------------------------------------------------------------------
  virtual LinearSolveStatus backsolve( const SparseVectorView_type& b, SparseVectorView_type& x ) const override;

protected:
  void init();

  using Base_type::A_;
  using Base_type::transpose_;
  AlgebraicEquationSetBase<DLA::MatrixD<Matrix_type>>& f_;
  std::vector< SparseMatrix_Diag<TM> > Dinv_;
};

} //namespace SLA
} //namespace SANS

#endif //ILU0_H
