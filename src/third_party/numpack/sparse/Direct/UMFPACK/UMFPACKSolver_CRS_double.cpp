// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define UMFPACK_INSTANTIATE
#include "UMFPACKSolver_impl.h"
#include "UMFPACKSolver_factorize.h"
#include "UMFPACKSolver_Solve_impl.h"

#include <memory>

namespace numpack
{
namespace SLA
{

//-----------------------------------------------------------------------------
template< class Matrix_type >
LinearSolveStatus UMFPACK<Matrix_type>::
backsolve( const int solvecode, const SparseVectorView_type& b, SparseVectorView_type& x ) const
{
  SANS_ASSERT( A_->m() == b.m() );
  SANS_ASSERT( A_->n() == x.m() );

  timer solvetime; // start timing solve

  SANS_UMFPACK_INT status = 0;

  status = SANS_UMFPACK_SOLVE(solvecode, Ap_, Ai_, Ax_, &x[0], &b[0], Numeric_, control_.data(), info_.data() );

  if ( status != UMFPACK_OK )
    throw( UMFPACKException(status, info_.data()) );

  if (timing_) std::cout << "UMFPACK solve time : " << solvetime.elapsed() << " second(s)" << std::endl;

  return LinearSolveStatus(true);
}

template class UMFPACK< SparseMatrix_CRS<double> >;

} //namespace SLA
} //namespace numpack
