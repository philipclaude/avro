// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef LINEARSOLVERBASE_H
#define LINEARSOLVERBASE_H

#include "tools/SANSException.h"

#include "numpack/VectorType.h"

namespace numpack 
{
namespace SLA
{

//----------------------------------------------------------------------------//
// A base class for all linear solvers and preconditioners
//----------------------------------------------------------------------------//

// All linear solvers solve the sparse system Ax = b, i.e. x = Inverse(A)*b
// May also solve the transposed system A'x = b

enum LinearSystemSolve
{
  RegularSolve,
  TransposeSolve
};

// Linear solve status returned
struct LinearSolveStatus
{
  explicit LinearSolveStatus(bool success) : success(success) {}
  LinearSolveStatus() : success(true) {}

  bool success;
};


template<class Matrix_type_>
class LinearSolverBase
{
public:
  typedef Matrix_type_ Matrix_type;
  typedef typename VectorType<Matrix_type>::type SparseVector_type;
  typedef typename VectorType<Matrix_type>::Viewtype SparseVectorView_type;

  LinearSolverBase(LinearSystemSolve solve) : A_(nullptr), transpose_(solve == TransposeSolve) {}
  virtual ~LinearSolverBase() {}

  // returns a reference to the matrix
  const Matrix_type& A() const { SANS_ASSERT(A_ != nullptr); return *A_; }

//-----------------------------------------------------------------------------
  // Factorizes the matrix given by an algebraic equation set
  virtual void factorize() = 0;

//-----------------------------------------------------------------------------
  // Solve the linear system (this assumes that the matrix has already been factorized)
  virtual LinearSolveStatus backsolve(const SparseVectorView_type& b, SparseVectorView_type& x) const = 0;

//-----------------------------------------------------------------------------
  virtual LinearSolveStatus solve(SparseVectorView_type& b, SparseVectorView_type& x)
  {
    // factorize the matrix and then solve
    factorize();
    return backsolve( b, x );
  }

  LinearSystemSolve systemSolve()
  {
    return transpose_ ? TransposeSolve : RegularSolve;
  }

protected:
  // the sparse matrix solved by this linear solver
  Matrix_type* A_;

  // Determines if solving the transposed system
  bool transpose_;
};

} //namespace SLA
} //namespace numpack 


#endif //LINEARSOLVERBASE_H
