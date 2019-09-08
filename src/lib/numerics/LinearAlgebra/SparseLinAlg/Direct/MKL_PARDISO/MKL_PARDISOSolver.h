// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MKL_PARDISOSOLVER_H
#define MKL_PARDISOSOLVER_H

#include "Python/PyDict.h"
#include "Python/Parameter.h"

#include "tools/SANSException.h"
#include "tools/noncopyable.h"

#include "LinearAlgebra/AlgebraicEquationSetBase.h"
#include "LinearAlgebra/SparseLinAlg/LinearSolverBase.h"
#include "LinearAlgebra/SparseLinAlg/SparseLinAlg_Inverse.h"
#include "LinearAlgebra/SparseLinAlg/ScalarMatrix_CRS.h"

#include <memory> // std::shared_ptr

#include "MKL_PARDISOSolver_defines.h"

namespace SANS
{
namespace SLA
{

//=============================================================================
struct MKL_PARDISOException : public SANSException
{
  explicit MKL_PARDISOException(const int status);

  virtual ~MKL_PARDISOException() throw() {}
};

//=============================================================================
//Forward declare
template< class Matrix_type >
class MKL_PARDISO;

//=============================================================================
struct MKL_PARDISOParam : noncopyable
{
  const ParameterBool Timing{"Timing", false, "Time Components of UMFPACK Solve"};

  template<class Matrix_type>
  static std::shared_ptr< LinearSolverBase<Matrix_type> >
  newSolver(const PyDict& SolverParam, AlgebraicEquationSetBase<Matrix_type>& f, LinearSystemSolve solve = RegularSolve)
  {
    typedef std::shared_ptr< LinearSolverBase<Matrix_type> > Solver_ptr;

    return Solver_ptr( new MKL_PARDISO<Matrix_type>( SolverParam, f, solve ) );
  }

  static void checkInputs(PyDict d);
  static MKL_PARDISOParam params;
};


//=============================================================================
template< class Matrix_type >
class MKL_PARDISO : public LinearSolverBase< Matrix_type >
{
public:
  typedef LinearSolverBase< Matrix_type > Base_type;

  typedef typename Base_type::SparseVector_type SparseVector_type;
  typedef typename Base_type::SparseVectorView_type SparseVectorView_type;
  typedef typename NonZeroPatternType<Matrix_type>::type SystemNonZeroPattern;

  MKL_PARDISOParam& params;

//-----------------------------------------------------------------------------
  MKL_PARDISO( AlgebraicEquationSetBase<Matrix_type>& f, LinearSystemSolve solve = RegularSolve, bool timing = false );
  MKL_PARDISO( const PyDict& d, AlgebraicEquationSetBase<Matrix_type>& f, LinearSystemSolve solve = RegularSolve );
  virtual ~MKL_PARDISO();

//-----------------------------------------------------------------------------
  virtual void factorize();

//-----------------------------------------------------------------------------
  virtual LinearSolveStatus backsolve( const SparseVectorView_type& b, SparseVectorView_type& x ) const;

//-----------------------------------------------------------------------------
  LinearSolveStatus backsolveTranspose( const SparseVectorView_type& b, SparseVectorView_type& x ) const;

protected:

//-----------------------------------------------------------------------------
  LinearSolveStatus backsolve( SANS_MKL_PARDISO_INT *iparm, const SparseVectorView_type& b, SparseVectorView_type& x ) const;

//-----------------------------------------------------------------------------
  void init();
  void deallocate();
  void free_pardiso();

protected:
  using Base_type::A_;
  using Base_type::transpose_;
  AlgebraicEquationSetBase<Matrix_type>& f_;
  ScalarMatrix_CRS<SANS_MKL_PARDISO_INT> *Ms_;         // Scalar CRS matrix passed to MKL PARDISO
  SANS_MKL_PARDISO_INT m_;                             // The number of rows in the matrix is needed to free the memory
  mutable SANS_MKL_PARDISO_INT mtype_;                 // Matrix type (symmetric, non-symmetric, etc.)
  mutable SANS_MKL_PARDISO_INT maxfct_;                // Maximum number of numerical factorizations.
  mutable SANS_MKL_PARDISO_INT mnum_;                  // Which factorization to use.
  mutable SANS_MKL_PARDISO_INT msglvl_;                // Output level
  mutable SANS_MKL_PARDISO_INT iparm_[64];             // Pardiso control parameters
  mutable void *pt_[64];                               // Internal solver memory pointer pt
  bool timing_;
};

} //namespace SLA
} //namespace SANS

#endif //MKL_PARDISOSOLVER_H
