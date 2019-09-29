// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "WritePlainVector.h"

namespace numpack 
{
namespace SLA
{

// write a plain (i.e. unformatted) vector to output stream
std::ostream& WritePlainVector( const ScalarVector& vec_plain, std::ostream& out )
{
  for (int j = 0; j < vec_plain.m; ++j)
    out << std::setprecision(16) << vec_plain.v[j] << std::endl;

  return out;
}

} // namespace SLA
} // namespace numpack 
