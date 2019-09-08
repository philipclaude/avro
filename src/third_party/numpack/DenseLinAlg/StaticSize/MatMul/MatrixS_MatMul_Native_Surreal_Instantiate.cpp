// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXS_MATMUL_NATIVE_INSTANTIATE
#include "MatrixS_MatMul_Native_impl.h"

#include "numpack/types/SurrealS.h"
#include "numpack/types/SurrealD.h"
#include "numpack/DenseLinAlg/StaticSize/VectorS.h"

namespace numpack 
{
namespace DLA
{
//-----------------------------------------------------------------------------
// Instantiates the expressions
//
// C = s * A x B
//
// where s is a scalar, and A, B, C are matrices. The template arguments are:
//
// MATRIXS_MATMUL_NATIVE( Am,An, Bm,Bn, A-type, B-type, s-type, C-type );
//
// Am, An are the sizes of the A matrix (similarly for B). The C matrix is deduced from A and B sizes.

//---- SurrealD = Real * SurrealD x SurrealD ------------------------------------
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealD, SurrealD, Real, SurrealD );
MATRIXS_MATMUL_NATIVE( 2,2, 2,2, SurrealD, SurrealD, Real, SurrealD );

//---- SurrealS = Real * Real x Real ----------------------------------------
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Real, Real, Real, SurrealS<1> );

MATRIXS_MATMUL_NATIVE( 2,1, 1,2, Real, Real, Real, SurrealS<2> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Real, Real, Real, SurrealS<2> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,2, Real, Real, Real, SurrealS<2> );

MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Real, Real, Real, SurrealS<3> );
MATRIXS_MATMUL_NATIVE( 3,3, 3,3, Real, Real, Real, SurrealS<3> );

MATRIXS_MATMUL_NATIVE( 4,4, 4,1, Real, Real, Real, SurrealS<4> );
MATRIXS_MATMUL_NATIVE( 4,4, 4,4, Real, Real, Real, SurrealS<4> );

MATRIXS_MATMUL_NATIVE( 5,5, 5,1, Real, Real, Real, SurrealS<5> );
MATRIXS_MATMUL_NATIVE( 5,5, 5,5, Real, Real, Real, SurrealS<5> );

MATRIXS_MATMUL_NATIVE( 2,2, 2,2, Real, Real, Real, SurrealS<6> );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Real, Real, Real, SurrealS<6> );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Real, Real, Real, SurrealS<7> );

MATRIXS_MATMUL_NATIVE( 6,6, 6,1, Real, Real, Real, SurrealS<6> );
MATRIXS_MATMUL_NATIVE( 7,7, 7,1, Real, Real, Real, SurrealS<7> );

//---- SurrealS = Real * Real x SurrealS ----------------------------------------
MATRIXS_MATMUL_NATIVE( 1,1, 1,1, Real, SurrealS<1>, Real, SurrealS<1> );
MATRIXS_MATMUL_NATIVE( 1,2, 2,2, Real, SurrealS<1>, Real, SurrealS<1> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Real, SurrealS<1>, Real, SurrealS<1> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,2, Real, SurrealS<1>, Real, SurrealS<1> );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Real, SurrealS<1>, Real, SurrealS<1> );

MATRIXS_MATMUL_NATIVE( 1,1, 1,1, Real, SurrealS<2>, Real, SurrealS<2> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Real, SurrealS<2>, Real, SurrealS<2> );

MATRIXS_MATMUL_NATIVE( 1,2, 2,2, Real, SurrealS<3>, Real, SurrealS<3> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Real, SurrealS<3>, Real, SurrealS<3> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,2, Real, SurrealS<3>, Real, SurrealS<3> );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Real, SurrealS<3>, Real, SurrealS<3> );
MATRIXS_MATMUL_NATIVE( 3,3, 3,3, Real, SurrealS<3>, Real, SurrealS<3> );

MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Real, SurrealS<4>, Real, SurrealS<4> );
MATRIXS_MATMUL_NATIVE( 4,4, 4,1, Real, SurrealS<4>, Real, SurrealS<4> );
MATRIXS_MATMUL_NATIVE( 4,4, 4,4, Real, SurrealS<4>, Real, SurrealS<4> );

MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Real, SurrealS<5>, Real, SurrealS<5> );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Real, SurrealS<5>, Real, SurrealS<5> );
MATRIXS_MATMUL_NATIVE( 5,5, 5,1, Real, SurrealS<5>, Real, SurrealS<5> );

MATRIXS_MATMUL_NATIVE( 1,2, 2,2, Real, SurrealS<6>, Real, SurrealS<6> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Real, SurrealS<6>, Real, SurrealS<6> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,2, Real, SurrealS<6>, Real, SurrealS<6> );
MATRIXS_MATMUL_NATIVE( 1,3, 3,3, Real, SurrealS<6>, Real, SurrealS<6> );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Real, SurrealS<6>, Real, SurrealS<6> );
MATRIXS_MATMUL_NATIVE( 3,3, 3,3, Real, SurrealS<6>, Real, SurrealS<6> );
MATRIXS_MATMUL_NATIVE( 6,6, 6,1, Real, SurrealS<6>, Real, SurrealS<6> );

MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Real, SurrealS<7>, Real, SurrealS<7> );
MATRIXS_MATMUL_NATIVE( 7,7, 7,1, Real, SurrealS<7>, Real, SurrealS<7> );

MATRIXS_MATMUL_NATIVE( 2,2, 2,2, Real, SurrealS<9>, Real, SurrealS<9> );

MATRIXS_MATMUL_NATIVE( 1,3, 3,3, Real, SurrealS<12>, Real, SurrealS<12> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,2, Real, SurrealS<12>, Real, SurrealS<12> );
MATRIXS_MATMUL_NATIVE( 3,3, 3,3, Real, SurrealS<12>, Real, SurrealS<12> );

MATRIXS_MATMUL_NATIVE( 1,4, 4,4, Real, SurrealS<20>, Real, SurrealS<20> );
MATRIXS_MATMUL_NATIVE( 4,4, 4,4, Real, SurrealS<20>, Real, SurrealS<20> );

MATRIXS_MATMUL_NATIVE( 3,3, 3,3, Real, SurrealS<24>, Real, SurrealS<24> );

MATRIXS_MATMUL_NATIVE( 3,3, 3,3, Real, SurrealS<48>, Real, SurrealS<48> );

typedef SurrealS<9,SurrealS<9>> SurrealS9_9;
typedef SurrealS<12,SurrealS<12>> SurrealS12_12;
typedef SurrealS<24,SurrealS<24>> SurrealS24_24;
typedef SurrealS<48,SurrealS<48>> SurrealS48_48;

MATRIXS_MATMUL_NATIVE( 2,2, 2,2, Real, SurrealS9_9  , Real, SurrealS9_9   );
MATRIXS_MATMUL_NATIVE( 2,2, 2,2, Real, SurrealS12_12, Real, SurrealS12_12 );

MATRIXS_MATMUL_NATIVE( 3,3, 3,3, Real, SurrealS24_24, Real, SurrealS24_24 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,3, Real, SurrealS48_48, Real, SurrealS48_48 );

//---- SurrealS = Real * int x SurrealS -----------------------------------------
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<1>, int, Real, SurrealS<1> );

//---- SurrealS = Real * SurrealS x Real ----------------------------------------
MATRIXS_MATMUL_NATIVE( 1,1, 1,1, SurrealS<1>, Real, Real, SurrealS<1> );
MATRIXS_MATMUL_NATIVE( 1,2, 2,1, SurrealS<1>, Real, Real, SurrealS<1> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<1>, Real, Real, SurrealS<1> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,2, SurrealS<1>, Real, Real, SurrealS<1> );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, SurrealS<1>, Real, Real, SurrealS<1> );
MATRIXS_MATMUL_NATIVE( 4,4, 4,1, SurrealS<1>, Real, Real, SurrealS<1> );

MATRIXS_MATMUL_NATIVE( 1,1, 1,1, SurrealS<2>, Real, Real, SurrealS<2> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<2>, Real, Real, SurrealS<2> );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, SurrealS<2>, Real, Real, SurrealS<2> );

MATRIXS_MATMUL_NATIVE( 1,1, 1,1, SurrealS<3>, Real, Real, SurrealS<3> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<3>, Real, Real, SurrealS<3> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,2, SurrealS<3>, Real, Real, SurrealS<3> );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, SurrealS<3>, Real, Real, SurrealS<3> );

MATRIXS_MATMUL_NATIVE( 1,1, 1,1, SurrealS<4>, Real, Real, SurrealS<4> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<4>, Real, Real, SurrealS<4> );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, SurrealS<4>, Real, Real, SurrealS<4> );
MATRIXS_MATMUL_NATIVE( 4,4, 4,1, SurrealS<4>, Real, Real, SurrealS<4> );

MATRIXS_MATMUL_NATIVE( 1,2, 2,2, SurrealS<5>, Real, Real, SurrealS<5> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<5>, Real, Real, SurrealS<5> );
MATRIXS_MATMUL_NATIVE( 5,5, 5,1, SurrealS<5>, Real, Real, SurrealS<5> );

MATRIXS_MATMUL_NATIVE( 1,2, 2,1, SurrealS<6>, Real, Real, SurrealS<6> );
MATRIXS_MATMUL_NATIVE( 1,2, 2,2, SurrealS<6>, Real, Real, SurrealS<6> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<6>, Real, Real, SurrealS<6> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,2, SurrealS<6>, Real, Real, SurrealS<6> );
MATRIXS_MATMUL_NATIVE( 1,3, 3,3, SurrealS<6>, Real, Real, SurrealS<6> );
MATRIXS_MATMUL_NATIVE( 3,3, 3,3, SurrealS<6>, Real, Real, SurrealS<6> );
MATRIXS_MATMUL_NATIVE( 6,6, 6,1, SurrealS<6>, Real, Real, SurrealS<6> );

MATRIXS_MATMUL_NATIVE( 1,3, 3,3, SurrealS<7>, Real, Real, SurrealS<7> );
MATRIXS_MATMUL_NATIVE( 7,7, 7,1, SurrealS<7>, Real, Real, SurrealS<7> );

MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<8>, Real, Real, SurrealS<8> );

MATRIXS_MATMUL_NATIVE( 1,3, 3,1, SurrealS<12>, Real, Real, SurrealS<12> );
MATRIXS_MATMUL_NATIVE( 3,3, 3,3, SurrealS<12>, Real, Real, SurrealS<12> );

MATRIXS_MATMUL_NATIVE( 4,4, 4,4, SurrealS<20>, Real, Real, SurrealS<20> );

//---- SurrealS = Real * SurrealS x SurrealS -----------------------------------
MATRIXS_MATMUL_NATIVE( 1,1, 1,1, SurrealS<1>, SurrealS<1>, Real, SurrealS<1> );
MATRIXS_MATMUL_NATIVE( 1,2, 2,1, SurrealS<1>, SurrealS<1>, Real, SurrealS<1> );
MATRIXS_MATMUL_NATIVE( 1,2, 2,2, SurrealS<1>, SurrealS<1>, Real, SurrealS<1> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<1>, SurrealS<1>, Real, SurrealS<1> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,2, SurrealS<1>, SurrealS<1>, Real, SurrealS<1> );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, SurrealS<1>, SurrealS<1>, Real, SurrealS<1> );
MATRIXS_MATMUL_NATIVE( 4,4, 4,1, SurrealS<1>, SurrealS<1>, Real, SurrealS<1> );
MATRIXS_MATMUL_NATIVE( 4,4, 4,2, SurrealS<1>, SurrealS<1>, Real, SurrealS<1> );
MATRIXS_MATMUL_NATIVE( 4,4, 4,3, SurrealS<1>, SurrealS<1>, Real, SurrealS<1> );
MATRIXS_MATMUL_NATIVE( 4,4, 4,4, SurrealS<1>, SurrealS<1>, Real, SurrealS<1> );
MATRIXS_MATMUL_NATIVE( 6,6, 6,1, SurrealS<1>, SurrealS<1>, Real, SurrealS<1> );

MATRIXS_MATMUL_NATIVE( 1,1, 1,1, SurrealS<2>, SurrealS<2>, Real, SurrealS<2> );
MATRIXS_MATMUL_NATIVE( 2,1, 1,2, SurrealS<2>, SurrealS<2>, Real, SurrealS<2> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<2>, SurrealS<2>, Real, SurrealS<2> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,2, SurrealS<2>, SurrealS<2>, Real, SurrealS<2> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,2, SurrealS<9>, SurrealS<9>, Real, SurrealS<9> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,2, SurrealS<12>, SurrealS<12>, Real, SurrealS<12> );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, SurrealS<2>, SurrealS<2>, Real, SurrealS<2> );

MATRIXS_MATMUL_NATIVE( 4,4, 4,1, SurrealS<2>, SurrealS<2>, Real, SurrealS<2> );
MATRIXS_MATMUL_NATIVE( 4,4, 4,2, SurrealS<2>, SurrealS<2>, Real, SurrealS<2> );
MATRIXS_MATMUL_NATIVE( 4,4, 4,3, SurrealS<2>, SurrealS<2>, Real, SurrealS<2> );
MATRIXS_MATMUL_NATIVE( 4,4, 4,4, SurrealS<2>, SurrealS<2>, Real, SurrealS<2> );

MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<3>, SurrealS<3>, Real, SurrealS<3> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,2, SurrealS<3>, SurrealS<3>, Real, SurrealS<3> );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, SurrealS<3>, SurrealS<3>, Real, SurrealS<3> );
MATRIXS_MATMUL_NATIVE( 3,3, 3,3, SurrealS<3>, SurrealS<3>, Real, SurrealS<3> );

MATRIXS_MATMUL_NATIVE( 4,4, 4,1, SurrealS<4>, SurrealS<4>, Real, SurrealS<4> );
MATRIXS_MATMUL_NATIVE( 4,4, 4,2, SurrealS<4>, SurrealS<4>, Real, SurrealS<4> );
MATRIXS_MATMUL_NATIVE( 4,4, 4,3, SurrealS<4>, SurrealS<4>, Real, SurrealS<4> );
MATRIXS_MATMUL_NATIVE( 4,4, 4,4, SurrealS<4>, SurrealS<4>, Real, SurrealS<4> );

MATRIXS_MATMUL_NATIVE( 5,5, 5,1, SurrealS<5>, SurrealS<5>, Real, SurrealS<5> );
MATRIXS_MATMUL_NATIVE( 5,5, 5,5, SurrealS<5>, SurrealS<5>, Real, SurrealS<5> );
MATRIXS_MATMUL_NATIVE( 5,5, 5,5, Real, SurrealS<5>, Real, SurrealS<5> );
MATRIXS_MATMUL_NATIVE( 5,5, 5,5, SurrealS<5>, Real, Real, SurrealS<5> );

MATRIXS_MATMUL_NATIVE( 2,1, 1,2, SurrealS<6>, SurrealS<6>, Real, SurrealS<6> );
MATRIXS_MATMUL_NATIVE( 1,2, 2,1, SurrealS<6>, SurrealS<6>, Real, SurrealS<6> );
MATRIXS_MATMUL_NATIVE( 1,2, 2,2, SurrealS<6>, SurrealS<6>, Real, SurrealS<6> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<6>, SurrealS<6>, Real, SurrealS<6> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,2, SurrealS<6>, SurrealS<6>, Real, SurrealS<6> );
MATRIXS_MATMUL_NATIVE( 3,3, 3,3, SurrealS<6>, SurrealS<6>, Real, SurrealS<6> );
MATRIXS_MATMUL_NATIVE( 6,6, 6,1, SurrealS<6>, SurrealS<6>, Real, SurrealS<6> );
MATRIXS_MATMUL_NATIVE( 6,6, 6,6, SurrealS<6>, SurrealS<6>, Real, SurrealS<6> );

MATRIXS_MATMUL_NATIVE( 7,7, 7,1, SurrealS<7>, SurrealS<7>, Real, SurrealS<7> );
MATRIXS_MATMUL_NATIVE( 7,7, 7,7, SurrealS<7>, SurrealS<7>, Real, SurrealS<7> );

MATRIXS_MATMUL_NATIVE( 2,1, 1,2, SurrealS<8>, SurrealS<8>, Real, SurrealS<8> );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<8>, SurrealS<8>, Real, SurrealS<8> );
MATRIXS_MATMUL_NATIVE( 8,8, 8,1, SurrealS<8>, SurrealS<8>, Real, SurrealS<8> );

MATRIXS_MATMUL_NATIVE( 1,2, 2,2, SurrealS<9>, SurrealS<9>, Real, SurrealS<9> );

MATRIXS_MATMUL_NATIVE( 5,5, 5,1, SurrealS<10>, SurrealS<10>, Real, SurrealS<10> );

MATRIXS_MATMUL_NATIVE( 3,3, 3,3, SurrealS<12>, SurrealS<12>, Real, SurrealS<12> );

MATRIXS_MATMUL_NATIVE( 4,4, 4,4, SurrealS<20>, SurrealS<20>, Real, SurrealS<20> );

MATRIXS_MATMUL_NATIVE( 3,3, 3,3, SurrealS<24>, SurrealS<24>, Real, SurrealS<24> );

MATRIXS_MATMUL_NATIVE( 3,3, 3,3, SurrealS<48>, SurrealS<48>, Real, SurrealS<48> );


MATRIXS_MATMUL_NATIVE( 2,2, 2,2, SurrealS9_9  , SurrealS9_9  , Real, SurrealS9_9   );
MATRIXS_MATMUL_NATIVE( 2,2, 2,2, SurrealS12_12, SurrealS12_12, Real, SurrealS12_12 );

MATRIXS_MATMUL_NATIVE( 3,3, 3,3, SurrealS24_24, SurrealS24_24, Real, SurrealS24_24 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,3, SurrealS48_48, SurrealS48_48, Real, SurrealS48_48 );


//---- SurrealS = SurrealS * Real x Real ------------------------------------
MATRIXS_MATMUL_NATIVE( 2,1, 1,2, Real, Real, SurrealS<2>, SurrealS<2> );
MATRIXS_MATMUL_NATIVE( 2,1, 1,2, Real, Real, SurrealS<3>, SurrealS<3> );
MATRIXS_MATMUL_NATIVE( 2,1, 1,2, Real, Real, SurrealS<4>, SurrealS<4> );
MATRIXS_MATMUL_NATIVE( 2,1, 1,2, Real, Real, SurrealS<5>, SurrealS<5> );
MATRIXS_MATMUL_NATIVE( 2,1, 1,2, Real, Real, SurrealS<6>, SurrealS<6> );

//---- SurrealS = SurrealS * SurrealS x SurrealS ------------------------------------
MATRIXS_MATMUL_NATIVE( 1,1, 1,1, SurrealS<1>, SurrealS<1>, SurrealS<1>, SurrealS<1> );
MATRIXS_MATMUL_NATIVE( 2,1, 1,2, SurrealS<1>, SurrealS<1>, SurrealS<1>, SurrealS<1> );
MATRIXS_MATMUL_NATIVE( 2,1, 1,2, SurrealS<3>, SurrealS<3>, SurrealS<3>, SurrealS<3> );
MATRIXS_MATMUL_NATIVE( 2,1, 1,2, SurrealS<5>, SurrealS<5>, SurrealS<5>, SurrealS<5> );
MATRIXS_MATMUL_NATIVE( 2,1, 1,2, SurrealS<6>, SurrealS<6>, SurrealS<6>, SurrealS<6> );
MATRIXS_MATMUL_NATIVE( 2,1, 1,2, SurrealS<8>, SurrealS<8>, SurrealS<8>, SurrealS<8> );

//-----------------------------------------------------------------------------------
// Multiplications of block matrices
//-----------------------------------------------------------------------------------
typedef VectorS<2,Real> Vector2;
typedef VectorS<2,SurrealS<2>> VectorS2;
typedef MatrixS<2,2,SurrealS<2>> MatrixS22;
MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS22, Vector2, Real, VectorS2 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS22, Vector2, Real, VectorS2 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS22, Vector2, Real, VectorS2 );

MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS22, VectorS2, Real, VectorS2 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS22, VectorS2, Real, VectorS2 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS22, VectorS2, Real, VectorS2 );

typedef VectorS<3,Real> Vector3;
typedef VectorS<3,SurrealS<3>> VectorS3;
typedef VectorS<3,SurrealS<1>> VectorS31;
typedef MatrixS<3,3,SurrealS<3>> MatrixS33;
typedef MatrixS<3,3,SurrealS<1>> MatrixS31;
MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS33, MatrixS33, Real, MatrixS33 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS33, MatrixS33, Real, MatrixS33 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS33, MatrixS33, Real, MatrixS33 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS33, MatrixS33, Real, MatrixS33 );

MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS33, VectorS3, Real, VectorS3 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS33, VectorS3, Real, VectorS3 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS33, VectorS3, Real, VectorS3 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS33, VectorS3, Real, VectorS3 );

MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS33, VectorS3, Real, VectorS3 );
MATRIXS_MATMUL_NATIVE( 1,3, 3,3, MatrixS33, VectorS3, Real, VectorS3 );
MATRIXS_MATMUL_NATIVE( 1,3, 3,1, MatrixS33, VectorS3, Real, VectorS3 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,3, MatrixS33, VectorS3, Real, VectorS3 );

MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS33, Vector3, Real, VectorS3 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS33, Vector3, Real, VectorS3 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS33, Vector3, Real, VectorS3 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS33, Vector3, Real, VectorS3 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS33, Vector3, Real, VectorS3 );
MATRIXS_MATMUL_NATIVE( 1,3, 3,3, MatrixS33, Vector3, Real, VectorS3 );
MATRIXS_MATMUL_NATIVE( 1,3, 3,1, MatrixS33, Vector3, Real, VectorS3 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,3, MatrixS33, Vector3, Real, VectorS3 );

MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS31, Vector3, Real, VectorS31 );

typedef VectorS<4,Real> Vector4;
typedef VectorS<4,SurrealS<4>> VectorS4;
typedef MatrixS<4,4,SurrealS<4>> MatrixS44;
MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS44, MatrixS44, Real, MatrixS44 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS44, MatrixS44, Real, MatrixS44 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS44, MatrixS44, Real, MatrixS44 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS44, MatrixS44, Real, MatrixS44 );

MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS44, VectorS4, Real, VectorS4 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS44, VectorS4, Real, VectorS4 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS44, VectorS4, Real, VectorS4 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS44, VectorS4, Real, VectorS4 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS44, VectorS4, Real, VectorS4 );

MATRIXS_MATMUL_NATIVE( 4,4, 4,1, MatrixS44, VectorS4, Real, VectorS4 );
MATRIXS_MATMUL_NATIVE( 1,4, 4,4, MatrixS44, VectorS4, Real, VectorS4 );
MATRIXS_MATMUL_NATIVE( 1,4, 4,1, MatrixS44, VectorS4, Real, VectorS4 );
MATRIXS_MATMUL_NATIVE( 4,4, 4,4, MatrixS44, VectorS4, Real, VectorS4 );

MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS44, Vector4, Real, VectorS4 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS44, Vector4, Real, VectorS4 );

typedef VectorS<5,Real> Vector5;
typedef VectorS<5,SurrealS<5>> VectorS5;
typedef MatrixS<5,5,SurrealS<5>> MatrixS55;

MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS55, Vector5, Real, VectorS5 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS55, Vector5, Real, VectorS5 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS55, Vector5, Real, VectorS5 );

MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS55, VectorS5, Real, VectorS5 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS55, VectorS5, Real, VectorS5 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS55, VectorS5, Real, VectorS5 );
MATRIXS_MATMUL_NATIVE( 5,5, 5,1, MatrixS55, VectorS5, Real, VectorS5 );

typedef VectorS<5,SurrealS<10>> VectorS10_5;
typedef MatrixS<5,5,SurrealS<10>> MatrixS10_55;

MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS10_55, Vector5, Real, VectorS10_5 );

MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS10_55, VectorS10_5, Real, VectorS10_5 );

typedef VectorS<6,Real> Vector6;
typedef VectorS<6,SurrealS<6>> VectorS6;
typedef MatrixS<6,6,SurrealS<6>> MatrixS66;

MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS66, Vector6, Real, VectorS6);
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS66, Vector6, Real, VectorS6 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS66, Vector6, Real, VectorS6 );

MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS66, VectorS6, Real, VectorS6 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS66, VectorS6, Real, VectorS6 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS66, VectorS6, Real, VectorS6 );

typedef VectorS<7,Real> Vector7;
typedef VectorS<7,SurrealS<7>> VectorS7;
typedef MatrixS<7,7,SurrealS<7>> MatrixS77;

MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS77, Vector7, Real, VectorS7 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS77, VectorS7, Real, VectorS7 );

typedef VectorS<8,SurrealS<8>> VectorS8;
typedef MatrixS<8,8,SurrealS<8>> MatrixS88;

MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS88, VectorS8, Real, VectorS8 );

// begin pcaplan
typedef MatrixS<4,4,Real> MatrixR44;
typedef MatrixS<4,4,SurrealS<10,Real>> MatrixS44_S10;
typedef MatrixS<1,4,Real> MatrixR14;
typedef MatrixS<1,4,SurrealS<10,Real>> MatrixS14_S10;
MATRIXS_MATMUL_NATIVE( 4,4,4,4, MatrixS44_S10, MatrixR44 , Real, MatrixS44_S10 );
MATRIXS_MATMUL_NATIVE( 4,4,4,4, MatrixS44_S10, MatrixS44_S10 , Real, MatrixS44_S10 );
MATRIXS_MATMUL_NATIVE( 4,4,4,4, MatrixR44, MatrixS44_S10 , Real, MatrixS44_S10 );
MATRIXS_MATMUL_NATIVE( 1,4,4,4, MatrixR14, MatrixS44_S10, Real, MatrixS14_S10 );
MATRIXS_MATMUL_NATIVE( 1,4,4,4, SurrealS<10>, SurrealS<10>, Real, SurrealS<10> );
MATRIXS_MATMUL_NATIVE( 4,4,4,4, SurrealS<10>, SurrealS<10>, Real, SurrealS<10> );
MATRIXS_MATMUL_NATIVE( 4,4,4,4, SurrealS<10>, Real, Real, SurrealS<10> );
MATRIXS_MATMUL_NATIVE( 4,4,4,4, Real, SurrealS<10>, Real, SurrealS<10> );
MATRIXS_MATMUL_NATIVE( 1,4,4,4, Real, SurrealS<10>, Real, SurrealS<10> );

MATRIXS_MATMUL_NATIVE( 4,4,4,4, Real,SurrealS<50>,Real,SurrealS<50>);
MATRIXS_MATMUL_NATIVE( 4,4,4,4, SurrealS<50>,SurrealS<50>,Real,SurrealS<50>);


// end pcaplan

} //namespace DLA
} //namespace numpack 
