// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SPARSELINALG_SCALAR_H
#define SPARSELINALG_SCALAR_H

#include "tools/SANSTraitsScalar.h"

namespace SANS
{

namespace SLA
{
  // Forward declare
  template<class T>
  class SparseVector;
}

template<class T>
struct Scalar< SLA::SparseVector<T> >
{
  typedef typename Scalar<T>::type type;
};

}

#endif //SPARSELINALG_SCALAR_H
