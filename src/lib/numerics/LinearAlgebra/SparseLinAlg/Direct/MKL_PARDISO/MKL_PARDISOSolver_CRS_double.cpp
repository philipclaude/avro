// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MKL_PARDISO_INSTANTIATE
#include <memory>
#include "MKL_PARDISOSolver_impl.h"
#include "MKL_PARDISOSolver_factorize.h"

namespace SANS
{
namespace SLA
{

//-----------------------------------------------------------------------------
template< class Matrix_type >
LinearSolveStatus MKL_PARDISO<Matrix_type>::backsolve( SANS_MKL_PARDISO_INT *iparm, const SparseVectorView_type& b, SparseVectorView_type& x ) const
{
  SANS_ASSERT( b.m() == x.m() );
  SANS_ASSERT( A_->m() == x.m() );

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
                   &nrhs, iparm, &msglvl_, (void*)&b[0], (void*)&x[0], &error);

  if ( error != 0 ) BOOST_THROW_EXCEPTION( MKL_PARDISOException(error) );

  if (timing_) std::cout << "MKL_PARDISO solve time : " << solvetime.elapsed() << " second(s)" << std::endl;

  return LinearSolveStatus(true);
}


template class MKL_PARDISO< SparseMatrix_CRS<double> >;

} //namespace SLA
} //namespace SANS
