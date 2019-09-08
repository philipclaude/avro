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
LinearSolveStatus UMFPACK<Matrix_type>::solve(SparseVectorView_type& b, SparseVectorView_type& x)
{
  if ( staticCondensed_ )
  {
    factorize(b, transpose_);

    SANS_ASSERT( x.m() == 3 && b.m() == 3 ); //hardcode to static condensation from 3x3 system
    SparseVectorView_type bcondensed = b.sub(1,2);
    SparseVectorView_type xcondensed = x.sub(1,2);

    LinearSolveStatus status = backsolve(bcondensed, xcondensed);

    f_.completeUpdate(b, xcondensed, x);

    return status;
  }
  else
  {
    factorize();
    return backsolve( b, x );
  }
}

} // namespace SLA
} // namespace numpack 
