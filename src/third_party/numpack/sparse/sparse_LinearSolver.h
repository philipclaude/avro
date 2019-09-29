// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef sparse_LINEARSOLVER_H
#define sparse_LINEARSOLVER_H

//Python must be included first
#include "Python/PyDict.h"

#include <memory>

#include "numpack/AlgebraicEquationSetBase.h"

#include "sparse_Type.h"
#include "LinearSolverBase.h"
#include "sparse_Inverse.h"

#include "Krylov/FGMRES.h"
#include "Direct/UMFPACK/UMFPACKSolver.h"
#ifdef SANS_PETSC
#include "PETSc/PETScSolver.h"
#endif
#ifdef INTEL_MKL
#include "Direct/MKL_PARDISO/MKL_PARDISOSolver.h"
#endif
namespace numpack 
{
namespace SLA
{

//=============================================================================
struct LinearSolverParam : noncopyable
{
  struct LinearSolverOptions
  {
    typedef DictKeyPair ExtractType;
    const ParameterString Solver{"Solver", "FGMRES", "Linear Solver Name" };
    const ParameterString& key = Solver;

    const DictOption FGMRES{"FGMRES", FGMRESParam::checkInputs};
    const DictOption UMFPACK{"UMFPACK", UMFPACKParam::checkInputs};
#ifdef SANS_PETSC
    const DictOption PETSc{"PETSc", PETScSolverParam::checkInputs};
#endif
#ifdef INTEL_MKL
    const DictOption MKL_PARDISO{"MKL_PARDISO", MKL_PARDISOParam::checkInputs};
#endif
    const std::vector<DictOption> options{ FGMRES
                                         , UMFPACK
#ifdef SANS_PETSC
                                         , PETSc
#endif
#ifdef INTEL_MKL
                                         , MKL_PARDISO
#endif
                                         };
  };
  const ParameterOption<LinearSolverOptions> LinearSolver{"LinearSolver", NO_DEFAULT, "Iterative or Direct Linear Solver"};

  static void checkInputs(PyDict d);
  static LinearSolverParam params;
};


//=============================================================================
template<class Matrix_type>
class LinearSolver
{
public:
  typedef std::shared_ptr< LinearSolverBase<Matrix_type> > Solver_ptr;
  typedef typename VectorType<Matrix_type>::Viewtype SparseVectorView_type;

  LinearSolverParam& params = LinearSolverParam::params;

  explicit LinearSolver( const PyDict& d, AlgebraicEquationSetBase<Matrix_type>& f, LinearSystemSolve solve = RegularSolve )
  {
    DictKeyPair SolverParam = d.get(params.LinearSolver);

    if ( SolverParam == params.LinearSolver.FGMRES )
      solver_ = FGMRESParam::newSolver<Matrix_type>( SolverParam, f, solve );
    else if ( SolverParam == params.LinearSolver.UMFPACK )
      solver_ = UMFPACKParam::newSolver<Matrix_type>( SolverParam, f, solve );
#ifdef SANS_PETSC
    else if ( SolverParam == params.LinearSolver.PETSc )
      solver_ = PETScSolverParam::newSolver<Matrix_type>( SolverParam, f, solve );
#endif
#ifdef INTEL_MKL
    else if ( SolverParam == params.LinearSolver.MKL_PARDISO )
      solver_ = MKL_PARDISOParam::newSolver<Matrix_type>( SolverParam, f, solve );
#endif
    else
      SANS_DEVELOPER_EXCEPTION("Unknown linear solver: %s", SolverParam.key().c_str());
  }

  const Matrix_type& A() { return solver_->A(); }

//-----------------------------------------------------------------------------
  void factorize()
  {
    //Factorize the matrix
    solver_->factorize();
  }

//-----------------------------------------------------------------------------
  LinearSolveStatus backsolve( const SparseVectorView_type& b, SparseVectorView_type& x) const
  {
    //This assumes that factorize has already been called
    return solver_->backsolve( b, x );
  }

//-----------------------------------------------------------------------------
  LinearSolveStatus
  solve( SparseVectorView_type& b, SparseVectorView_type& x )
  {
    //This invokes both the factorization and the back solve
    return solver_->solve( b, x );
  }

  void dumpJacobian(std::string filename)
  {
    std::cout << "Writing Jacobian matrix to file: " << filename << "..." << std::endl;
    WriteMatrixMarketFile( solver_->A(), filename );
  }

private:
  Solver_ptr solver_;
};


} //namespace SLA
} //namespace numpack 

#endif //sparse_LINEARSOLVER_H
