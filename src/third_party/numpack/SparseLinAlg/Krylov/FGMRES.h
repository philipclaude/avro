// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef FGMRES_H
#define FGMRES_H

#include "Python/PyDict.h" //Python must be first
#include "Python/Parameter.h"

#include "numpack/SparseLinAlg/LinearSolverBase.h"
#include "numpack/SparseLinAlg/Preconditioners/Identity.h"
#include "numpack/SparseLinAlg/Preconditioners/LU_SGS.h"

#include <memory>

namespace numpack 
{
namespace SLA
{
//Forward declare
template< class Matrix_type >
class FGMRES;

//=============================================================================
struct FGMRESParam : noncopyable
{
  const ParameterNumeric<int> nOuter{"nOuter", 20 , 0, NO_LIMIT, "Number of Outer Iterations"};
  const ParameterNumeric<int> nInner{"nInner", 100, 0, NO_LIMIT, "Number of Inner Iterations"};

  const ParameterNumeric<Real> tol{"tol", 1e-12, 0, NO_LIMIT, "Convergence Tolerance"};

  const ParameterBool PrintCovergence{"PrintConvergence", false, "Print Convergence History"};

  struct PreconditionerOptions
  {
    typedef DictKeyPair ExtractType;
    const ParameterString Name{"Name", NO_DEFAULT, "Name of the preconditioner" };
    const ParameterString& key = Name;

    const DictOption Identity{"Identity", IdentityParam::checkInputs};
    const DictOption LU_SGS  {"LU_SGS"  , LU_SGSParam::checkInputs  };

    const std::vector<DictOption> options{Identity, LU_SGS};
  };
  const ParameterOption<PreconditionerOptions> Preconditioner{"Preconditioner", NO_DEFAULT, "Preconditioner for FGMRES"};


  template<class Matrix_type>
  static std::shared_ptr< LinearSolverBase<Matrix_type> >
  newSolver( const PyDict& SolverParam, AlgebraicEquationSetBase<Matrix_type>& f, LinearSystemSolve solve = RegularSolve )
  {
    typedef std::shared_ptr< LinearSolverBase<Matrix_type> > Solver_ptr;

    DictKeyPair PrecondParam = SolverParam.get(params.Preconditioner);

    if ( PrecondParam == params.Preconditioner.Identity )
      return Solver_ptr( new FGMRES<Matrix_type >( SolverParam, Solver_ptr( new Identity<Matrix_type>(f, solve) ) ) );

    if ( PrecondParam == params.Preconditioner.LU_SGS )
      return Solver_ptr( new FGMRES<Matrix_type >( SolverParam, Solver_ptr( new LU_SGS<Matrix_type>(f, solve) ) ) );

    SANS_ASSERT( false );
    return Solver_ptr(); //Just to avoid compiler warning
  }

  static void checkInputs(PyDict d);
  static FGMRESParam params;
};


//=============================================================================
template< class Matrix_type >
class FGMRES : public LinearSolverBase< Matrix_type >
{
public:
  typedef LinearSolverBase< Matrix_type > Base_type;
  typedef std::shared_ptr< LinearSolverBase<Matrix_type> > Solver_ptr;

  typedef typename Base_type::SparseVector_type SparseVector_type;
  typedef typename Base_type::SparseVectorView_type SparseVectorView_type;

  FGMRESParam& params;

protected:
  //Norm_type norm;

  using Base_type::A_; //The sparse matrix
  Solver_ptr M_;       //The preconditioner

  Real tol_;               // Tolerance for convergence
  mutable Real tolFinal_;  // Final tolerance
  int nInner_;             // Maximum number of inner iterations
  int nOuter_;             // Maximum number of outer iterations
  bool printConv_;         // Print convergence history
  mutable int itTotal_;

public:
//-----------------------------------------------------------------------------
  FGMRES( const PyDict& d, Solver_ptr M );
  virtual ~FGMRES();

//-----------------------------------------------------------------------------
  virtual void factorize() override;

//-----------------------------------------------------------------------------
  virtual LinearSolveStatus backsolve( const SparseVectorView_type& b, SparseVectorView_type& x ) const override;

  using Base_type::backsolve;
};

} //namespace SLA
} //namespace numpack 

#endif //FGMRES_H
