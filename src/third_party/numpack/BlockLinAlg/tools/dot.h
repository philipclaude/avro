// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SPARSELINALG_DOT_H
#define SPARSELINALG_DOT_H

#include "tools/SANSnumerics.h"
#include "numpack/SparseLinAlg/SparseLinAlg_Type.h"
#include "SparseLinAlg_Scalar.h"

#include "numpack/DenseLinAlg/tools/norm.h"

namespace numpack 
{

template<class T>>
inline Real
dot(const BLA::VectorBlock_2<DLA::VectorD<SLA::SparseVector<T> >, DLA::VectorD<SLA::SparseVector<T> > >& eL,
    const BLA::VectorBlock_2<DLA::VectorD<SLA::SparseVector<T> >, DLA::VectorD<SLA::SparseVector<T> > >& eR)
{
  SANS_ASSERT(eL.v0.m() == eR.v0.m());
  SANS_ASSERT(eL.v1.m() == eR.v1.m());

  return dot(eL.v0, eR.v0) + dot(eL.v1, eR.v1);
}


}

#endif //SPARSELINALG_DOT_H
