// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXD_DECOMPOSE_LU_INSTANTIATE
#include "MatrixD_Decompose_LU_impl.h"

namespace tinymat 
{
namespace DLA
{

typedef MatrixS<1,1,Real> MatrixS11;
typedef MatrixS<2,2,Real> MatrixS22;
typedef MatrixS<3,3,Real> MatrixS33;
typedef MatrixS<4,4,Real> MatrixS44;
typedef MatrixS<5,5,Real> MatrixS55;
typedef MatrixS<6,6,Real> MatrixS66;

template struct MatrixDLU< MatrixS11 >;
template struct MatrixDLU< MatrixS22 >;
template struct MatrixDLU< MatrixS33 >;
template struct MatrixDLU< MatrixS44 >;
template struct MatrixDLU< MatrixS55 >;
template struct MatrixDLU< MatrixS66 >;

template struct MatrixDLU< MatrixS<2,2,MatrixS<2,2,Real>> >;
template struct MatrixDLU< MatrixS<2,2,MatrixS<3,3,Real>> >;
template struct MatrixDLU< MatrixS<2,2,MatrixS<4,4,Real>> >;
template struct MatrixDLU< MatrixS<2,2,MatrixS<5,5,Real>> >;

template struct MatrixDLU< MatrixS<3,3,MatrixS<2,2,Real>> >;
template struct MatrixDLU< MatrixS<3,3,MatrixS<3,3,Real>> >;
template struct MatrixDLU< MatrixS<3,3,MatrixS<4,4,Real>> >;
template struct MatrixDLU< MatrixS<3,3,MatrixS<5,5,Real>> >;
template struct MatrixDLU< MatrixS<3,3,MatrixS<6,6,Real>> >;

}
}
