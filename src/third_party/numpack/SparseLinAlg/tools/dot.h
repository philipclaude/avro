// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SPARSELINALG_DOT_H
#define SPARSELINALG_DOT_H

#include "tools/SANSnumerics.h"
#include "numpack/SparseLinAlg/SparseLinAlg_Type.h"

namespace numpack 
{

template<class ExprL, class ExprR>
inline Real
dot(const SLA::SparseLinAlgType<ExprL, false>& eL, const SLA::SparseLinAlgType<ExprR, false>& eR)
{
  const ExprL& a = eL.cast();
  const ExprR& b = eR.cast();

  SANS_ASSERT(a.m() == b.m());

  const unsigned int size = a.m();

  Real val = 0;
  for (unsigned int i = 0; i < size; i++)
    val += dot(a[i],b[i]);

  return val;
}


}

#endif //SPARSELINALG_DOT_H
