// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "UpperTriangMatrix.h"

#include <iostream>
#include <iomanip>

namespace SANS
{
namespace SLA
{

//=============================================================================
std::ostream& operator<<(std::ostream& out, UpperTriangMatrix const &TriM) // output
{
  unsigned int i,j;
  for (i = 0; i < TriM.n_; i++)
  {
    for (j = 0; j < i; j++)
      out << std::setw(15) << "";
    for (j = i; j < TriM.n_; j++)
      out << std::setw(15) << std::setprecision(6) << TriM(i,j);
    out << std::endl;
  }
  return out;
}

} // namespace SLA
} // namespace SANS
