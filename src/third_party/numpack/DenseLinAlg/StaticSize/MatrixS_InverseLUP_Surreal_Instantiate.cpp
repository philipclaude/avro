// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXS_INVERSELUP_INSTANTIATE
#include "MatrixS_InverseLUP_impl.h"

#include "numpack/types/SurrealS.h"

namespace numpack 
{
namespace DLA
{

// Explicitly instantiate all datatypes used here

// A is Surreal, b is Surreal
template struct MatrixSLUPSolver< 2, 2, SurrealS<6,Real>, MatrixS<2,2,SurrealS<6,Real>> >;

template struct MatrixSLUPSolver< 3, 3, SurrealS<1,Real>, MatrixS<3,3,SurrealS<1,Real>> >;
template struct MatrixSLUPSolver< 3, 3, SurrealS<3,Real>, MatrixS<3,1,SurrealS<3,Real>> >;
template struct MatrixSLUPSolver< 3, 3, SurrealS<3,Real>, MatrixS<3,3,SurrealS<3,Real>> >;

template struct MatrixSLUPSolver< 4, 4, SurrealS<4,Real>, MatrixS<4,1,SurrealS<4,Real>> >;
template struct MatrixSLUPSolver< 4, 4, SurrealS<4,Real>, MatrixS<4,4,SurrealS<4,Real>> >;

template struct MatrixSLUPSolver< 5, 5, SurrealS<5,Real>, MatrixS<5,1,SurrealS<5,Real>> >;
template struct MatrixSLUPSolver< 5, 5, SurrealS<5,Real>, MatrixS<5,5,SurrealS<5,Real>> >;

template struct MatrixSLUPSolver< 6, 6, SurrealS<6,Real>, MatrixS<6,1,SurrealS<6,Real>> >;
template struct MatrixSLUPSolver< 6, 6, SurrealS<6,Real>, MatrixS<6,6,SurrealS<6,Real>> >;

template struct MatrixSLUPSolver< 7, 7, SurrealS<7,Real>, MatrixS<7,1,SurrealS<7,Real>> >;
template struct MatrixSLUPSolver< 7, 7, SurrealS<7,Real>, MatrixS<7,7,SurrealS<7,Real>> >;

// A is Real, b is Surreal
template struct MatrixSLUPSolver< 3, 3, Real, MatrixS<3,1,SurrealS<3,Real>> >;
template struct MatrixSLUPSolver< 4, 4, Real, MatrixS<4,1,SurrealS<4,Real>> >;

template struct MatrixSLUPSolver< 5, 5, Real, MatrixS<5,1,SurrealS<5,Real>> >;

template struct MatrixSLUPSolver< 6, 6, Real, MatrixS<6,1,SurrealS<6,Real>> >;

template struct MatrixSLUPSolver< 7, 7, Real, MatrixS<7,1,SurrealS<7,Real>> >;

}
}
