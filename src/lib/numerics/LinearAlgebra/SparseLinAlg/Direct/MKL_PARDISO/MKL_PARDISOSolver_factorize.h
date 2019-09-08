// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(MKL_PARDISO_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include <mkl_pardiso.h>

#include "MKL_PARDISOSolver.h"

#include "tools/timer.h"

namespace SANS
{
namespace SLA
{
//-----------------------------------------------------------------------------
template< class Matrix_type >
void MKL_PARDISO<Matrix_type>::init()
{
  if (timing_) std::cout << "Non-zero time : " << std::flush;

  timer nztime; // start timing non-zero pattern evaluation

  SystemNonZeroPattern nz(f_.matrixSize());
  f_.jacobian(nz);

  if (timing_) std::cout << nztime.elapsed() << " second(s)" << std::endl;

  // construct the matrix
  A_ = new Matrix_type(nz);
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
void MKL_PARDISO<Matrix_type>::factorize()
{
  if (timing_) std::cout << "Jacobian time : " << std::flush;
  timer jactime; // start timing Jacobian evaluation

  // compute the new matrix
  *A_ = 0;
  f_.jacobian(*A_);

  if (timing_) std::cout << jactime.elapsed() << " second(s)" << std::endl;

  delete Ms_;
  Ms_ = new ScalarMatrix_CRS<SANS_MKL_PARDISO_INT>(*A_);

  //Only solving square matricies
  SANS_ASSERT( Ms_->m() == Ms_->n() );
  m_ = Ms_->m();

  SANS_MKL_PARDISO_INT perm = 0;
  SANS_MKL_PARDISO_INT nrhs = 0;
  SANS_MKL_PARDISO_INT error = 0;
  SANS_MKL_PARDISO_INT phase = 0;
  double b = 0, x = 0;

  if (timing_) std::cout << "Numeric factor time : " << std::flush;
  timer numerictime;

  /* ----------------------------------------------------------*/
  /*   Reordering with Symbolic and Numeric Factorization.     */
  /*    This step also allocates all memory that is necessary  */
  /*    for the factorization.                                 */
  /* ----------------------------------------------------------*/
  phase = 12;
  SANS_MKL_PARDISO(pt_, &maxfct_, &mnum_, &mtype_, &phase,
                   &m_, Ms_->Rx(), Ms_->Rp(), Ms_->Ri(), &perm,
                   &nrhs, iparm_, &msglvl_, &b, &x, &error);

  if ( error != 0 ) BOOST_THROW_EXCEPTION( MKL_PARDISOException(error) );

  if (timing_) std::cout << numerictime.elapsed() << " second(s)" << std::endl;
}

} //namespace SLA
} //namespace SANS
