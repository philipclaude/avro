// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXS_DECOMPOSE_LU_INSTANTIATE
#include "MatrixS_Decompose_LU_impl.h"

#include "tinymat/types/SurrealS.h"

namespace tinymat 
{
namespace DLA
{

//Explicit instantiation
template struct MatrixSLU<2,SurrealS<6,Real>>;

template struct MatrixSLU<3,SurrealS<1,Real>>;
template struct MatrixSLU<4,SurrealS<4,Real>>;
template struct MatrixSLU<5,SurrealS<5,Real>>;
template struct MatrixSLU<6,SurrealS<6,Real>>;
template struct MatrixSLU<7,SurrealS<7,Real>>;
}
}
