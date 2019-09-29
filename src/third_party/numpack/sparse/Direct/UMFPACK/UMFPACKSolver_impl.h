// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(UMFPACK_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "UMFPACKSolver.h"

namespace numpack 
{
namespace SLA
{
//-----------------------------------------------------------------------------
template< class Matrix_type >
UMFPACK<Matrix_type>::
UMFPACK( AlgebraicEquationSetBase<Matrix_type>& f, LinearSystemSolve solve, bool timing ) :
  Base_type(solve),
  params(UMFPACKParam::params),
  f_(f),
  Symbolic_(nullptr),
  Numeric_(nullptr),
  Ap_(nullptr),
  Ai_(nullptr),
  Ax_(nullptr),
  timing_(timing),
  staticCondensed_(f.isStaticCondensed())
{
  // Set UMFPACK "sys" parameters here
  if (transpose_)
    solvecode_ = UMFPACK_At; // solve transposed system
  else
    solvecode_ = UMFPACK_A;  // solve regular system

  // Set UMFPACK "Control" parameter here
  SANS_UMFPACK_DEFAULTS(control_.data() ); // use UMFPACK default settings

  // "relative pivot tolerance for threshold partial pivoting with row interchanges"
  //   1.0 for TRUE partial pivoting (most accurate but expensive);
  //   <=0.0 for NO partial pivoting at all (least accurate but cheaper)
  control_[UMFPACK_PIVOT_TOLERANCE] = 1.0;

  // Clear UMFPACK "Control"; this is a bit redundant because umfpack_*_symbolic will do this at its start
  for (double& param : info_) { param = -1.0; }

  // Other initialization
  init();
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
UMFPACK<Matrix_type>::
UMFPACK( const PyDict& d,
         AlgebraicEquationSetBase<Matrix_type>& f,
         LinearSystemSolve solve ) :
  UMFPACK(f, solve, d.get(UMFPACKParam::params.Timing))
  {}

//-----------------------------------------------------------------------------
template< class Matrix_type >
LinearSolveStatus UMFPACK<Matrix_type>::
backsolve( const SparseVectorView_type& b, SparseVectorView_type& x ) const
{
  return backsolve(solvecode_, b, x);
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
LinearSolveStatus UMFPACK<Matrix_type>::
backsolveTranspose( const SparseVectorView_type& b, SparseVectorView_type& x ) const
{
  return backsolve(UMFPACK_At, b, x);
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
UMFPACK<Matrix_type>::~UMFPACK()
{
  delete A_;  A_  = nullptr;
  delete [] Ap_; Ap_ = nullptr;
  delete [] Ai_; Ai_ = nullptr;
  delete [] Ax_; Ax_ = nullptr;

  // deallocate UMFPACK's Numeric and Symbolic objects
  SANS_UMFPACK_FREE_NUMERIC(&Numeric_);
  SANS_UMFPACK_FREE_SYMBOLIC(&Symbolic_);
}

} //namespace SLA
} //namespace numpack 
