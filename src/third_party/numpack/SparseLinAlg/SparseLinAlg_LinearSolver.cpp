// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "SparseLinAlg_LinearSolver.h"

#define PYDICT_INSTANTIATE
#include "Python/PyDict_impl.h"

namespace numpack 
{

PYDICT_PARAMETER_OPTION_INSTANTIATE(SLA::LinearSolverParam::LinearSolverOptions)

namespace SLA
{

// cppcheck-suppress passedByValue
void LinearSolverParam::checkInputs(PyDict d)
{
  std::vector<const ParameterBase*> allParams;
  allParams.push_back(d.checkInputs(params.LinearSolver));
  d.checkUnknownInputs(allParams);
}
LinearSolverParam LinearSolverParam::params;


} //namespace SLA
} //namespace numpack 
