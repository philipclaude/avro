// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MKL_PARDISO_INSTANTIATE
#include "numpack/dense/dynamic/MatrixD.h"
#include "MKL_PARDISOSolver_factorize.h"
#include "MKL_PARDISOSolver_impl.h"

namespace numpack 
{
namespace SLA
{

//-----------------------------------------------------------------------------
template< class Matrix_type >
LinearSolveStatus MKL_PARDISO<Matrix_type>::backsolve( SANS_MKL_PARDISO_INT *iparm, const SparseVectorView_type& b, SparseVectorView_type& x ) const
{
  SANS_ASSERT( b.m() == x.m() );
  SANS_ASSERT( A_->m() == x.m() );

#ifdef INTEL_MKL

  timer solvetime; // start timing solve

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
                   &nrhs, iparm, &msglvl_, (void*)&b[0][0], (void*)&x[0][0], &error);

  if ( error != 0 ) BOOST_THROW_EXCEPTION( MKL_PARDISOException(error) );

  if (timing_) std::cout << "MKL_PARDISO solve time : " << solvetime.elapsed() << " second(s)" << std::endl;

#else

  BOOST_THROW_EXCEPTION( MKL_PARDISOException(42) );

#endif

  return LinearSolveStatus(true);
}

// TODO: Non-zero pattern fix needed

//template class MKL_PARDISO< SparseMatrix_CRS< DLA::MatrixD<Real> > >;

} //namespace SLA
} //namespace numpack 
