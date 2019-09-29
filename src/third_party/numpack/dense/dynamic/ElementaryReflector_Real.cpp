// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define ELEMENTARYREFLECTOR_INSTANTIATE
#include "ElementaryReflector_impl.h"

#include "numpack/dense/static/VectorS.h"

namespace numpack 
{
namespace DLA
{

template struct ElementaryReflector_impl<Real,Real>;

template struct ElementaryReflector_impl<Real, VectorS<1,Real> >;

}
}
