// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef sparse_DOT_H
#define sparse_DOT_H

#include "tools/SANSnumerics.h"
#include "numpack/sparse/sparse_Type.h"

namespace numpack 
{

template<class ExprL, class ExprR>
inline Real
dot(const SLA::sparseType<ExprL, false>& eL, const SLA::sparseType<ExprR, false>& eR)
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

#endif //sparse_DOT_H
