// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef UMFPACKSOLVER_H
#define UMFPACKSOLVER_H

#include <array> // std::array
#include <memory> // std::shared_ptr

#include "Python/PyDict.h"
#include "Python/Parameter.h"

#include "tools/SANSException.h"
#include "tools/noncopyable.h"

#include "numpack/AlgebraicEquationSetBase.h"

#include "numpack/sparse/LinearSolverBase.h"

#include "UMFPACKSolver_defines.h"

namespace numpack 
{
namespace SLA
{

//=============================================================================
struct UMFPACKException : public SANSException
{
  explicit UMFPACKException(const int status);
  UMFPACKException(const int status, const double *info);

  virtual ~UMFPACKException() throw() {}
};

//=============================================================================
//Forward declare
template< class Matrix_type >
class UMFPACK;

//=============================================================================
struct UMFPACKParam : noncopyable
{
  const ParameterBool Timing{"Timing", false, "Time Components of UMFPACK Solve"};

  template<class Matrix_type>
  static std::shared_ptr< LinearSolverBase<Matrix_type> >
  newSolver(const PyDict& SolverParam, AlgebraicEquationSetBase<Matrix_type>& f, LinearSystemSolve solve = RegularSolve )
  {
    typedef std::shared_ptr< LinearSolverBase<Matrix_type> > Solver_ptr;

    return Solver_ptr( new UMFPACK<Matrix_type>( SolverParam, f, solve ) );
  }

  static void checkInputs(PyDict d);
  static UMFPACKParam params;
};


//=============================================================================
template< class Matrix_type >
class UMFPACK : public LinearSolverBase< Matrix_type >
{
public:
  typedef LinearSolverBase< Matrix_type > Base_type;

  typedef typename Base_type::SparseVector_type SparseVector_type;
  typedef typename Base_type::SparseVectorView_type SparseVectorView_type;
  typedef typename NonZeroPatternType<Matrix_type>::type SystemNonZeroPattern;

  typedef std::array<double, UMFPACK_CONTROL> ControlParamType;
  typedef std::array<double, UMFPACK_INFO> InfoParamType;

  UMFPACKParam& params;

//-----------------------------------------------------------------------------
  UMFPACK( AlgebraicEquationSetBase<Matrix_type>& f, LinearSystemSolve solve = RegularSolve, bool timing = false );
  UMFPACK( const PyDict& d, AlgebraicEquationSetBase<Matrix_type>& f, LinearSystemSolve solve = RegularSolve );
  virtual ~UMFPACK();

//-----------------------------------------------------------------------------
  virtual void factorize() override;

  //-----------------------------------------------------------------------------
  void factorize( SparseVectorView_type& bcondensed, bool transpose); //for static condensation

//-----------------------------------------------------------------------------
  virtual LinearSolveStatus backsolve( const SparseVectorView_type& b, SparseVectorView_type& x ) const override;

//-----------------------------------------------------------------------------
  LinearSolveStatus backsolveTranspose( const SparseVectorView_type& b, SparseVectorView_type& x ) const;

  //-----------------------------------------------------------------------------
  virtual LinearSolveStatus solve(SparseVectorView_type& b, SparseVectorView_type& x) override;

  using Base_type::backsolve;

protected:
//-----------------------------------------------------------------------------
  void init();
  void factorizeMatrix();

//-----------------------------------------------------------------------------
  LinearSolveStatus backsolve( const int solvecode, const SparseVectorView_type& b, SparseVectorView_type& x ) const;

  using Base_type::A_;
  using Base_type::transpose_;
  AlgebraicEquationSetBase<Matrix_type>& f_;
  void *Symbolic_;
  void *Numeric_;
  // (Ap_, Ai_, Ax_) = (column_pointer, row_index, nonzero values)
  // forms CSC (compressed sparse column) sparse matrix storage which is used in UMFPACK
  SANS_UMFPACK_INT *Ap_, *Ai_;
  double *Ax_;
  int solvecode_; // UMFPACK "sys" parameter which determines which type of linear system is to be solved. (first argument to SANS_UMFPACK_SOLVE)
  ControlParamType control_; // store UMFPACK "Control" parameter.
  mutable InfoParamType info_;  // store UMFPACK "Info" parameter.
  bool timing_;
  bool staticCondensed_;
};

} //namespace SLA
} //namespace numpack 

#endif //UMFPACKSOLVER_H
