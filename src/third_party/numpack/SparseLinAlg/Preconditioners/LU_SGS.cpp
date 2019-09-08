// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "LU_SGS.h"


namespace numpack 
{
namespace SLA
{

// cppcheck-suppress passedByValue
void LU_SGSParam::checkInputs(PyDict d)
{
  std::vector<const ParameterBase*> allParams;
  d.checkUnknownInputs(allParams);
}
LU_SGSParam LU_SGSParam::params;


} //namespace SLA
} //namespace numpack 
