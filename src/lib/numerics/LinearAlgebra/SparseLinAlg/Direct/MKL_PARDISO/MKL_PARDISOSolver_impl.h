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



namespace SANS
{
namespace SLA
{

//-----------------------------------------------------------------------------
template< class Matrix_type >
MKL_PARDISO<Matrix_type>::MKL_PARDISO( AlgebraicEquationSetBase<Matrix_type>& f, LinearSystemSolve solve, bool timing ) :
  Base_type(solve),
  params(MKL_PARDISOParam::params),
  f_(f),
  Ms_(NULL),
  m_(0),
  mtype_(11),           // Real and nonsymmetric
  maxfct_(1),           // Maximum number of numerical factorizations.
  mnum_(1),             // Which factorization to use.
  msglvl_(0),           // No output
  timing_(timing)
{
  for (int i = 0; i < 64; i++ )
  {
    iparm_[i] = 0;
    pt_[i] = 0;
  }

  //It is recommended to use iparm[10] = 1 (scaling) and iparm[12]= 1
  //(matching) for highly indefinite symmetric matrices, for example
  //from interior point optimizations or saddle point problems.

  iparm_[0] = 1;         // No solver default
  iparm_[1] = 2;         // Fill-in reordering from METIS
  iparm_[3] = 0;         // No iterative-direct algorithm
  iparm_[4] = 0;         // No user fill-in reducing permutation
  iparm_[5] = 0;         // Write solution into x
  iparm_[7] = 3;         // Max numbers of iterative refinement steps
  iparm_[9] = 13;        // Perturb the pivot elements with 1E-13
  iparm_[10] = 1;        // Use nonsymmetric permutation and scaling MPS
  iparm_[11] = 0;        // Solve with transposed or conjugate transposed matrix A
  iparm_[12] = 1;        // Maximum weighted matching algorithm
  iparm_[13] = 0;        // Output: Number of perturbed pivots
  iparm_[17] = 0;        // Output: Number of nonzeros in the factor LU
  iparm_[18] = 0;        // Output: Mflops for LU factorization
  iparm_[19] = 0;        // Output: Numbers of CG Iterations
  iparm_[34] = 1;        // PARDISO use C-style indexing for ia and ja arrays

  if (transpose_)
    iparm_[11] = 2; // solve transposed system

  init();
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
MKL_PARDISO<Matrix_type>::MKL_PARDISO( const PyDict& d, AlgebraicEquationSetBase<Matrix_type>& f, LinearSystemSolve solve ) :
  MKL_PARDISO(f, solve, d.get(MKL_PARDISOParam::params.Timing))
{
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
LinearSolveStatus MKL_PARDISO<Matrix_type>::
backsolve( const SparseVectorView_type& b, SparseVectorView_type& x ) const
{
  return backsolve(iparm_, b, x);
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
LinearSolveStatus MKL_PARDISO<Matrix_type>::
backsolveTranspose( const SparseVectorView_type& b, SparseVectorView_type& x ) const
{
  SANS_MKL_PARDISO_INT iparm[64];
  for (int i = 0; i < 64; i++ )
    iparm[i] = iparm_[i];

  iparm[11] = 2; // solve transposed system

  return backsolve(iparm, b, x);
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
MKL_PARDISO<Matrix_type>::~MKL_PARDISO()
{
  deallocate();
  free_pardiso();
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
void MKL_PARDISO<Matrix_type>::deallocate()
{
  delete A_; A_ = NULL;
  delete Ms_; Ms_ = NULL;
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
void MKL_PARDISO<Matrix_type>::free_pardiso()
{
#ifdef INTEL_MKL
  SANS_MKL_PARDISO_INT error = 0;
  SANS_MKL_PARDISO_INT idum = 0; // Dummy
  double ddum = 0;    // Dummy

  //Release the internal memory
  SANS_MKL_PARDISO_INT phase = -1;
  SANS_MKL_PARDISO(pt_, &maxfct_, &mnum_, &mtype_, &phase,
                   &m_, &ddum, &idum, &idum, &idum, &idum,
                   iparm_, &msglvl_, &ddum, &ddum, &error);

  if ( error != 0 ) BOOST_THROW_EXCEPTION( MKL_PARDISOException(error) );

  //Reset pointers
  for (int i = 0; i < 64; i++ )
    pt_[i] = 0;

#else

  BOOST_THROW_EXCEPTION( MKL_PARDISOException(42) );

#endif
}


} //namespace SLA
} //namespace SANS
