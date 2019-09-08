// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(UMFPACK_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "UMFPACKSolver.h"
#include "LinearAlgebra/SparseLinAlg/ScalarVector.h"

#include "tools/timer.h"

namespace SANS
{
namespace SLA
{

//-----------------------------------------------------------------------------
template< class Matrix_type >
LinearSolveStatus UMFPACK<Matrix_type>::
backsolve( const int solvecode, const SparseVectorView_type& b, SparseVectorView_type& x ) const
{
  SANS_ASSERT( b.m() == x.m() );
  SANS_ASSERT( A_->m() == x.m() );

  ScalarVector xx(x);
  ScalarVector bb(b);

  timer solvetime; // start timing solve

  SANS_ASSERT( bb.m == xx.m );

  int status = SANS_UMFPACK_SOLVE(solvecode, Ap_, Ai_, Ax_, xx, bb, Numeric_, control_.data(), info_.data() );

  if ( status != UMFPACK_OK )
    BOOST_THROW_EXCEPTION( UMFPACKException(status, info_.data()) );

  if (timing_) std::cout << "UMFPACK solve time : " << solvetime.elapsed() << " second(s)" << std::endl;

#if 0
  //TODO: add an option to double check that the linear solve is actually successful
  //      because cases exist where no failure is reported but the linear solution is wrong!
#endif
  xx.setTo(x);

  return LinearSolveStatus(true);
}


}
}
