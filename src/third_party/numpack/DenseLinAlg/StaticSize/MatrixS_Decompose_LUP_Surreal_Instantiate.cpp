// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXS_DECOMPOSE_LUP_INSTANTIATE
#include "MatrixS_Decompose_LUP_impl.h"

#include "numpack/types/SurrealS.h"

namespace numpack 
{
namespace DLA
{

//Explicit instantiation
template struct MatrixSLUP<2,SurrealS<6,Real>>;

template struct MatrixSLUP<3,SurrealS<1,Real>>;
template struct MatrixSLUP<3,SurrealS<3,Real>>;
template struct MatrixSLUP<4,SurrealS<4,Real>>;
template struct MatrixSLUP<5,SurrealS<5,Real>>;
template struct MatrixSLUP<6,SurrealS<6,Real>>;
template struct MatrixSLUP<7,SurrealS<7,Real>>;
}
}
