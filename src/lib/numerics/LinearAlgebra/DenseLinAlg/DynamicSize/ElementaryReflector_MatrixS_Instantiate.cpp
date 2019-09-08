// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define ELEMENTARYREFLECTOR_INSTANTIATE
#include "ElementaryReflector_impl.h"

#include "LinearAlgebra/DenseLinAlg/StaticSize/MatrixS.h"

namespace SANS
{
namespace DLA
{

template struct ElementaryReflector_impl< MatrixS<1,1,Real>, VectorS<1,Real> >;
template struct ElementaryReflector_impl< MatrixS<2,2,Real>, VectorS<2,Real> >;
template struct ElementaryReflector_impl< MatrixS<4,4,Real>, VectorS<4,Real> >;

template struct ElementaryReflector_impl< MatrixS<1,1,Real>, MatrixS<1,1,Real> >;
template struct ElementaryReflector_impl< MatrixS<2,2,Real>, MatrixS<2,2,Real> >;
template struct ElementaryReflector_impl< MatrixS<4,4,Real>, MatrixS<4,4,Real> >;

}
}
