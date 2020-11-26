// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "SingularException.h"

#include <sstream>
#include <iomanip>

namespace tinymat 
{
namespace DLA
{

void
SingularMatrixException::errMessage(const Real denom)
{
  errString += "Singular Matrix!\n\n";
  std::stringstream msg;
  msg << std::scientific << std::setprecision(16);
  msg << "denominator = " << denom;
  errString += msg.str();
}


} // namespace DLA
} // namespace tinymat 
