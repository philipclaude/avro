// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "Identity.h"


namespace SANS
{
namespace SLA
{

// cppcheck-suppress passedByValue
void IdentityParam::checkInputs(PyDict d)
{
  std::vector<const ParameterBase*> allParams;
  d.checkUnknownInputs(allParams);
}
IdentityParam IdentityParam::params;


} //namespace SLA
} //namespace SANS
