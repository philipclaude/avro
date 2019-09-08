// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(MKL_PARDISO_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#ifdef INTEL_MKL
#include <mkl_pardiso.h>
#endif

#include "MKL_PARDISOSolver.h"
#include "LinearAlgebra/SparseLinAlg/ScalarVector.h"

#include "tools/timer.h"

namespace SANS
{
namespace SLA
{

//-----------------------------------------------------------------------------
template< class Matrix_type >
LinearSolveStatus MKL_PARDISO<Matrix_type>::backsolve( SANS_MKL_PARDISO_INT *iparm, const SparseVectorView_type& b, SparseVectorView_type& x ) const
{

#ifdef INTEL_MKL

  SANS_ASSERT( b.m() == x.m() );
  SANS_ASSERT( A_->m() == x.m() );

  ScalarVector xx(x);
  ScalarVector bb(b);

  timer solvetime; // start timing solve

  SANS_ASSERT( bb.m == xx.m );

  SANS_MKL_PARDISO_INT perm = 0;
  SANS_MKL_PARDISO_INT nrhs = 1;
  SANS_MKL_PARDISO_INT error = 0;
  SANS_MKL_PARDISO_INT phase = 0;

  SANS_MKL_PARDISO_INT m = Ms_->m();

  /* -----------------------------------------------*/
  /*  Back substitution and iterative refinement   */
  /* -----------------------------------------------*/

  phase = 33;
  SANS_MKL_PARDISO(pt_, &maxfct_, &mnum_, &mtype_, &phase,
                   &m, Ms_->Rx(), Ms_->Rp(), Ms_->Ri(), &perm,
                   &nrhs, iparm, &msglvl_, bb, xx, &error);

  if ( error != 0 ) BOOST_THROW_EXCEPTION( MKL_PARDISOException(error) );

  if (timing_) std::cout << "MKL_PARDISO solve time : " << solvetime.elapsed() << " second(s)" << std::endl;

  xx.setTo(x);

#else

  BOOST_THROW_EXCEPTION( MKL_PARDISOException(42) );

#endif

  return LinearSolveStatus(true);
}


}
}
