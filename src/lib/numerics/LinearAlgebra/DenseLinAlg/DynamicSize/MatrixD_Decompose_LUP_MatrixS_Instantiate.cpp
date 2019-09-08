// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXD_DECOMPOSE_LUP_INSTANTIATE
#include "MatrixD_Decompose_LUP_impl.h"

namespace SANS
{
namespace DLA
{

template struct MatrixDLUP< MatrixS<1,1,Real> >;
template struct MatrixDLUP< MatrixS<2,2,Real> >;
template struct MatrixDLUP< MatrixS<4,4,Real> >;
template struct MatrixDLUP< MatrixS<5,5,Real> >;
template struct MatrixDLUP< MatrixS<6,6,Real> >;
template struct MatrixDLUP< MatrixS<7,7,Real> >;

}
}
