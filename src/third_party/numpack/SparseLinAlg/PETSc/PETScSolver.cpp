// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "PETScSolver.h"

#include <petscoptions.h>
#include <petscerror.h>

#define PYDICT_INSTANTIATE
#include "Python/PyDict_impl.h"

namespace numpack 
{
PYDICT_PARAMETER_OPTION_INSTANTIATE(SLA::PreconditionerCommonParam::PreconditionerSideOptions)
PYDICT_PARAMETER_OPTION_INSTANTIATE(SLA::SubPreconditionerParam::PreconditionerOptions)
PYDICT_PARAMETER_OPTION_INSTANTIATE(SLA::SubPreconditionerParam::PreconditionerSideOptions)

PYDICT_PARAMETER_OPTION_INSTANTIATE(SLA::PETScSolverParam::PreconditionerOptions)
PYDICT_PARAMETER_OPTION_INSTANTIATE(SLA::PETScSolverParam::KSPSolverOptions)
PYDICT_PARAMETER_OPTION_INSTANTIATE(SLA::PreconditionerILUParam::OrderingOptions)

namespace SLA
{

//=============================================================================
PETScSolverException::PETScSolverException( const PetscErrorCode status )
{
  errString += "PETSc Error\n";

  const char *error_message = nullptr;
  PetscErrorMessage(status, &error_message, NULL);

  if (error_message == nullptr)
    errString += "No PETSc error message available";
  else
    errString += error_message;
}

void PreconditionerILUParam::checkInputs(PyDict d)
{
  std::vector<const ParameterBase*> allParams;
  allParams.push_back(d.checkInputs(params.FillLevel));
  allParams.push_back(d.checkInputs(params.NonzeroDiagonal));
  allParams.push_back(d.checkInputs(params.Ordering));
  allParams.push_back(d.checkInputs(params.PreconditionerSide));
  allParams.push_back(d.checkInputs(params.MDFOuterBlockSize));
  d.checkUnknownInputs(allParams);
}
PreconditionerILUParam PreconditionerILUParam::params;


void PreconditionerBlockJacobiParam::checkInputs(PyDict d)
{
  std::vector<const ParameterBase*> allParams;
  allParams.push_back(d.checkInputs(params.SubPreconditioner));
  allParams.push_back(d.checkInputs(params.PreconditionerSide));
  allParams.push_back(d.checkInputs(params.Ordering));
  d.checkUnknownInputs(allParams);
}
PreconditionerBlockJacobiParam PreconditionerBlockJacobiParam::params;

void PreconditionerASMParam::checkInputs(PyDict d)
{
  std::vector<const ParameterBase*> allParams;
  allParams.push_back(d.checkInputs(params.Overlap));
  allParams.push_back(d.checkInputs(params.SubPreconditioner));
  allParams.push_back(d.checkInputs(params.PreconditionerSide));
  allParams.push_back(d.checkInputs(params.Ordering));
  d.checkUnknownInputs(allParams);
}
PreconditionerASMParam PreconditionerASMParam::params;

void PreconditionerJacobiParam::checkInputs(PyDict d)
{
  std::vector<const ParameterBase*> allParams;
  allParams.push_back(d.checkInputs(params.PreconditionerSide));
  allParams.push_back(d.checkInputs(params.Ordering));
  d.checkUnknownInputs(allParams);
}
PreconditionerJacobiParam PreconditionerJacobiParam::params;

void PreconditionerHYPREParam::checkInputs(PyDict d)
{
  std::vector<const ParameterBase*> allParams;
  allParams.push_back(d.checkInputs(params.PreconditionerSide));
  allParams.push_back(d.checkInputs(params.Ordering));
  d.checkUnknownInputs(allParams);
}
PreconditionerHYPREParam PreconditionerHYPREParam::params;

void PreconditionerDirectLUParam::checkInputs(PyDict d)
{
  std::vector<const ParameterBase*> allParams;
  allParams.push_back(d.checkInputs(params.PreconditionerSide));
  allParams.push_back(d.checkInputs(params.Ordering));
  d.checkUnknownInputs(allParams);
}
PreconditionerDirectLUParam PreconditionerDirectLUParam::params;

void PreconditionerNoneParam::checkInputs(PyDict d)
{
  std::vector<const ParameterBase*> allParams;
  allParams.push_back(d.checkInputs(params.PreconditionerSide));
  allParams.push_back(d.checkInputs(params.Ordering));
  d.checkUnknownInputs(allParams);
}
PreconditionerNoneParam PreconditionerNoneParam::params;

void PETScSolverParam::checkInputs(PyDict d)
{
  PyDict copy(d);

  // check all other parameters
  std::vector<const ParameterBase*> allParams;
  allParams.push_back(d.checkInputs(params.RelativeTolerance));
  allParams.push_back(d.checkInputs(params.AbsoluteTolerance));
  allParams.push_back(d.checkInputs(params.DivergenceTolerance));
  allParams.push_back(d.checkInputs(params.MaxIterations));
  allParams.push_back(d.checkInputs(params.GMRES_Restart));
  allParams.push_back(d.checkInputs(params.Preconditioner));
  allParams.push_back(d.checkInputs(params.Verbose));
  allParams.push_back(d.checkInputs(params.Timing));
  allParams.push_back(d.checkInputs(params.Memory));
  allParams.push_back(d.checkInputs(params.computeSingularValues));
  allParams.push_back(d.checkInputs(params.printMatrixInfo));
  allParams.push_back(d.checkInputs(params.ResidualHistoryFile));
  allParams.push_back(d.checkInputs(params.FilenameBase));
  allParams.push_back(d.checkInputs(params.PETScOptions));
  allParams.push_back(d.checkInputs(params.PETScOptionsFile));
  allParams.push_back(d.checkInputs(params.PETScMatPrefix));
  allParams.push_back(d.checkInputs(params.PETScKSPPrefix));
  allParams.push_back(d.checkInputs(params.KSPSolver));
  d.checkUnknownInputs(allParams);

  std::string PETScOptions = copy.get(params.PETScOptions);
  std::string PETScOptionsFile = copy.get(params.PETScOptionsFile);

  // check that all PETSc options can be inserted
  PetscOptions options;
  PETSc_STATUS( PetscOptionsCreate(&options) );
  PETSc_STATUS( PetscOptionsInsertString(options, PETScOptions.c_str()) );
  if (!PETScOptionsFile.empty())
  {
    //TODO: Maybe something other than PETSC_COMM_SELF here...
    PETSc_STATUS( PetscOptionsInsertFile(PETSC_COMM_SELF, options, PETScOptionsFile.c_str(), PETSC_TRUE) );
  }
  PETSc_STATUS( PetscOptionsDestroy(&options) );
}
PETScSolverParam PETScSolverParam::params;

} //namespace SLA
} //namespace numpack 
