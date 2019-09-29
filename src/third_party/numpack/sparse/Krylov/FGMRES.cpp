// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "FGMRES.h"

#define PYDICT_INSTANTIATE
#include "Python/PyDict_impl.h"

namespace numpack 
{

PYDICT_PARAMETER_OPTION_INSTANTIATE(SLA::FGMRESParam::PreconditionerOptions)

namespace SLA
{

// cppcheck-suppress passedByValue
void FGMRESParam::checkInputs(PyDict d)
{
  std::vector<const ParameterBase*> allParams;
  allParams.push_back(d.checkInputs(params.nOuter));
  allParams.push_back(d.checkInputs(params.nInner));

  allParams.push_back(d.checkInputs(params.tol));

  allParams.push_back(d.checkInputs(params.PrintCovergence));

  allParams.push_back(d.checkInputs(params.Preconditioner));

  d.checkUnknownInputs(allParams);
}
FGMRESParam FGMRESParam::params;

} //namespace SLA
} //namespace numpack 
