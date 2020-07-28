// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(UMFPACK_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "UMFPACKSolver.h"
#include "numpack/sparse/ScalarMatrix_CRS.h"

#include <algorithm> // std::copy
#include <umfpack.h>

#undef USE_FACTORIZE_VERBOSE // turn this on to display verbose info from UMFPACK

#include "tools/timer.h"

namespace numpack
{
namespace SLA
{

//-----------------------------------------------------------------------------
template< class Matrix_type >
void UMFPACK<Matrix_type>::init()
{
  SANS_UMFPACK_INT status = 0;
  SANS_UMFPACK_INT m = 0;

  {
    if (timing_) std::cout << "Non-zero time : " << std::flush;

    timer nztime; // start timing non-zero pattern evaluation

    SystemNonZeroPattern nz(f_.matrixSize());

    if (staticCondensed_)
    {
      //dummy residual vector
      SparseVector_type dummyb(f_.vectorEqSize());
      dummyb = 0;

      f_.jacobian(dummyb, nz, transpose_);
    }
    else
    {
      f_.jacobian(nz);
    }

    if (timing_) std::cout << nztime.elapsed() << " second(s)" << std::endl;

    // construct the matrix
    A_ = new Matrix_type(nz);

    // clear the matrix
    (*A_) = 0;
  }

  {
#ifdef USE_FACTORIZE_VERBOSE
    timer scalartime;
#endif

    ScalarMatrix_CRS<SANS_UMFPACK_INT> Ms(*A_);

#ifdef USE_FACTORIZE_VERBOSE
    std::cout << "UMFPACK<>::factorize: scalar conversion time: " << scalartime.elapsed() << std::endl;
    std::cout << "UMFPACK<>::factorize: Ms.m = " << Ms.m() << "  Ms.nnz() = " << Ms.nnz() << std::endl;
#endif

    SANS_ASSERT( Ms.m() == Ms.n() ); // Can only solve square matrices

    m = Ms.m();

    Ap_ = new SANS_UMFPACK_INT[m+1];
    Ai_ = new SANS_UMFPACK_INT[Ms.nnz()];
    Ax_ = new double[Ms.nnz()];

    //Call transpose to convert the storage format to CSC
    status = SANS_UMFPACK_TRANSPOSE(m, m, Ms.Rp(), Ms.Ri(), Ms.Rx(), nullptr, nullptr, Ap_, Ai_, Ax_);

#ifdef USE_FACTORIZE_VERBOSE
    std::cout << "UMFPACK<>::factorize: status[SANS_UMFPACK_TRANSPOSE] = " << status << std::endl;
#endif

    if ( status != UMFPACK_OK )
      throw(UMFPACKException(status));
  }

  // perform symbolic factorization
  {
    if (timing_) std::cout << "Symbolic factor time : " << std::flush;

    timer symbolictime;

    // Ax_ is "Used only for gathering statistics about how many nonzeros are placed on the diagonal by the fill-reducing ordering"
    // but Ax_ is 0 because we have not yet computed the Jacobian, so it is not included here
    status = SANS_UMFPACK_SYMBOLIC(m, m, Ap_, Ai_, nullptr, &Symbolic_, control_.data(), info_.data());

#ifdef USE_FACTORIZE_VERBOSE
  std::cout << "UMFPACK<>::factorize: status[SANS_UMFPACK_SYMBOLIC] = " << status << std::endl;
#endif

    if ( status != UMFPACK_OK )
      throw( UMFPACKException(status, info_.data()) );

    if (timing_) std::cout << symbolictime.elapsed() << " second(s)" << std::endl;
  }

}

//-----------------------------------------------------------------------------
template< class Matrix_type >
void UMFPACK<Matrix_type>::
factorize(SparseVectorView_type& b, bool transpose)
{
  SANS_ASSERT(staticCondensed_ == true);
  timer jactime; // start timing Jacobian evaluation

  // update the matrix
  *A_ = 0;

//  const std::vector<std::vector<Real>> nrmRsd1 = f_.residualNorm(b);
//
//  std::string file1 = "tmp/rsd1_rank";
//  file1 += std::to_string(continuousmap_[0].comm->rank());
//  file1 += ".dat";
//  WritePlainVector( b, file1 );

  f_.jacobian(b, *A_, transpose);

//  const std::vector<std::vector<Real>> nrmRsd2 = f_.residualNorm(b);
//
//  std::string file2 = "tmp/rsd2_rank";
//  file2 += std::to_string(continuousmap_[0].comm->rank());
//  file2 += ".dat";
//  WritePlainVector( b, file2 );

  if (timing_) std::cout << jactime.elapsed() << " second(s)" << std::endl;

  factorizeMatrix();
}


//-----------------------------------------------------------------------------
template< class Matrix_type >
void UMFPACK<Matrix_type>::factorize()
{
  SANS_ASSERT(staticCondensed_ == false);
  timer jactime; // start timing Jacobian evaluation

  // update the matrix
  *A_ = 0;
  f_.jacobian(*A_);

  if (timing_) std::cout << jactime.elapsed() << " second(s)" << std::endl;

  factorizeMatrix();
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
void UMFPACK<Matrix_type>::factorizeMatrix()
{
  SANS_UMFPACK_INT status = 0;

  {
    // create the scalar matrix (TODO: This is a HACK)
    ScalarMatrix_CRS<SANS_UMFPACK_INT> Ms(*A_);

#if 0
    std::string filename = "tmp/jac.mtx";
    std::cout << "Writing Jacobian matrix to file: " << filename << "..." << std::endl;
    Ms.WriteMatrixMarketFile( filename );
#endif

    //Call transpose to convert the sparse matrix storage format to CSC
    status = SANS_UMFPACK_TRANSPOSE(Ms.m(), Ms.m(), Ms.Rp(), Ms.Ri(), Ms.Rx(), nullptr, nullptr, Ap_, Ai_, Ax_);

    if ( status != UMFPACK_OK )
      throw( UMFPACKException(status) );
  }

  timer numerictime;
  if (timing_) std::cout << "Numeric factor time : " << std::flush;

  if (Numeric_ != nullptr) { SANS_UMFPACK_FREE_NUMERIC(&Numeric_); } // Deallocate any lingering Numeric object
  status = SANS_UMFPACK_NUMERIC(Ap_, Ai_, Ax_, Symbolic_, &Numeric_, control_.data(), info_.data());

#ifdef USE_FACTORIZE_VERBOSE
  std::cout << "UMFPACK<>::factorize: numeric time: " << numerictime.elapsed() << std::endl;
  std::cout << "UMFPACK<>::factorize: status[SANS_UMFPACK_NUMERIC] = " << status << std::endl;

  ControlParamType control_verbose;
  std::copy(control_.begin(), control_.end(), control_verbose.begin() );
  control_verbose[UMFPACK_PRL] = 2;

  std::cout << "--------------- Display SANS_UMFPACK_REPORT_INFO ---------------" << std::endl;
  SANS_UMFPACK_REPORT_INFO( control_verbose.data(), info_.data() );
#endif

  if ( status != UMFPACK_OK )
    throw( UMFPACKException(status, info_.data()) );

  if (timing_) std::cout << numerictime.elapsed() << " second(s)" << std::endl;
}

} //namespace SLA
} //namespace numpack
