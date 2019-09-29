// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXD_DECOMPOSE_LU_INSTANTIATE
#include "MatrixD_Decompose_LU_impl.h"

namespace numpack 
{
namespace DLA
{

template struct MatrixDLU< MatrixD<Real> >;

template struct MatrixDLU< MatrixD<MatrixS<1,1,Real>> >;
template struct MatrixDLU< MatrixD<MatrixS<2,2,Real>> >;
template struct MatrixDLU< MatrixD<MatrixS<3,3,Real>> >;
template struct MatrixDLU< MatrixD<MatrixS<4,4,Real>> >;
template struct MatrixDLU< MatrixD<MatrixS<5,5,Real>> >;
template struct MatrixDLU< MatrixD<MatrixS<6,6,Real>> >;
template struct MatrixDLU< MatrixD<MatrixS<7,7,Real>> >;

}
}
