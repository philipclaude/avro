// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXD_MATMUL_NATIVE_INSTANTIATE
#include "MatMul_Native_impl.h"

// Allow instantiation of MatrixS matmul here to improve performance
#define MATRIXS_MATMUL_NATIVE_INSTANTIATE
#include "numpack/DenseLinAlg/StaticSize/MatMul/MatrixS_MatMul_Native_impl.h"

namespace numpack 
{
namespace DLA
{
//-----------------------------------------------------------------------------
// Instantiates the expressions
//
// C = A x B
//
// where s is a scalar, and A, B, C are matrices. The template arguments are:
//
// MatMul_Native< A-type, B-type, C-type >
//
template class MatMul_Native< Real, MatrixS<1,8,Real>, MatrixS<1,8,Real> >;

template class MatMul_Native< MatrixS<1,1,Real>, MatrixS<1,1,Real>, MatrixS<1,1,Real> >;
template class MatMul_Native< MatrixS<2,2,Real>, MatrixS<2,2,Real>, MatrixS<2,2,Real> >;
template class MatMul_Native< MatrixS<3,3,Real>, MatrixS<3,3,Real>, MatrixS<3,3,Real> >;
template class MatMul_Native< MatrixS<4,4,Real>, MatrixS<4,4,Real>, MatrixS<4,4,Real> >;
template class MatMul_Native< MatrixS<5,5,Real>, MatrixS<5,5,Real>, MatrixS<5,5,Real> >;
template class MatMul_Native< MatrixS<6,6,Real>, MatrixS<6,6,Real>, MatrixS<6,6,Real> >;
template class MatMul_Native< MatrixS<7,7,Real>, MatrixS<7,7,Real>, MatrixS<7,7,Real> >;

template class MatMul_Native< MatrixS<1,8,Real>, MatrixS<8,1,Real>, Real >;
template class MatMul_Native< MatrixS<2,8,Real>, MatrixS<8,1,Real>, MatrixS<2,1,Real> >; // e.g. coupled IBL2/panel without transition/lag
template class MatMul_Native< MatrixS<4,8,Real>, MatrixS<8,1,Real>, MatrixS<4,1,Real> >; // e.g. coupled IBL2/panel with transition/lag

template class MatMul_Native< MatrixS<2,1,Real>, MatrixS<1,2,Real>, MatrixS<2,2,Real> >;
template class MatMul_Native< MatrixS<3,1,Real>, MatrixS<1,3,Real>, MatrixS<3,3,Real> >;
template class MatMul_Native< MatrixS<4,1,Real>, MatrixS<1,4,Real>, MatrixS<4,4,Real> >;
template class MatMul_Native< MatrixS<5,1,Real>, MatrixS<1,5,Real>, MatrixS<5,5,Real> >;
template class MatMul_Native< MatrixS<6,1,Real>, MatrixS<1,6,Real>, MatrixS<6,6,Real> >;
template class MatMul_Native< MatrixS<7,1,Real>, MatrixS<1,7,Real>, MatrixS<7,7,Real> >;

template class MatMul_Native< MatrixS<8,4,Real>, MatrixS<4,1,Real>, MatrixS<8,1,Real> >; // e.g. coupled IBL2/panel with transition/lag
template class MatMul_Native< MatrixS<8,8,Real>, MatrixS<8,1,Real>, MatrixS<8,1,Real> >;

template class MatMul_Native< MatrixS<1,1,Real>, VectorS<1,Real>, VectorS<1,Real> >;
template class MatMul_Native< MatrixS<2,2,Real>, VectorS<2,Real>, VectorS<2,Real> >;
template class MatMul_Native< MatrixS<3,3,Real>, VectorS<3,Real>, VectorS<3,Real> >;
template class MatMul_Native< MatrixS<4,4,Real>, VectorS<4,Real>, VectorS<4,Real> >;
template class MatMul_Native< MatrixS<5,5,Real>, VectorS<5,Real>, VectorS<5,Real> >;
template class MatMul_Native< MatrixS<6,6,Real>, VectorS<6,Real>, VectorS<6,Real> >;
template class MatMul_Native< MatrixS<7,7,Real>, VectorS<7,Real>, VectorS<7,Real> >;

template class MatMul_Native< MatrixS<1,1,Real>, VectorS<1,Real>, Real >;
template class MatMul_Native< MatrixS<1,2,Real>, VectorS<2,Real>, Real >;
template class MatMul_Native< MatrixS<1,3,Real>, VectorS<3,Real>, Real >;
template class MatMul_Native< MatrixS<1,4,Real>, VectorS<4,Real>, Real >;

typedef VectorS<2,Real> VectorS2;
typedef VectorS<3,Real> VectorS3;
typedef VectorS<4,Real> VectorS4;
typedef VectorS<5,Real> VectorS5;
typedef VectorS<6,Real> VectorS6;
typedef VectorS<7,Real> VectorS7;

typedef MatrixS<2,2,Real> MatrixS22;
typedef MatrixS<3,2,Real> MatrixS32;
typedef MatrixS<3,3,Real> MatrixS33;
typedef MatrixS<4,4,Real> MatrixS44;
typedef MatrixS<5,5,Real> MatrixS55;
typedef MatrixS<6,6,Real> MatrixS66;
typedef MatrixS<7,7,Real> MatrixS77;

template class MatMul_Native< MatrixS<1,1,Real>, MatrixS<1,1,VectorS3>, MatrixS<1,1,VectorS3> >;
template class MatMul_Native< MatrixS<1,1,VectorS3>, MatrixS<1,1,Real>, MatrixS<1,1,VectorS3> >;

template class MatMul_Native< MatrixS<1,1,MatrixS33>, VectorS<1,MatrixS33>, MatrixS33 >;
template class MatMul_Native< MatrixS<1,2,MatrixS33>, VectorS<2,MatrixS33>, MatrixS33 >;

// -- DGBR2 mtxPDEElemL_rL*rL_qL etc.
template class MatMul_Native< MatrixS<1,1,MatrixS22>, VectorS<1,Real>, MatrixS22 >;
template class MatMul_Native< MatrixS<1,2,MatrixS22>, VectorS<2,Real>, MatrixS22 >;
template class MatMul_Native< MatrixS<1,3,MatrixS22>, VectorS<3,Real>, MatrixS22 >;
template class MatMul_Native< MatrixS<1,1,MatrixS33>, VectorS<1,Real>, MatrixS33 >;
template class MatMul_Native< MatrixS<1,2,MatrixS33>, VectorS<2,Real>, MatrixS33 >;
template class MatMul_Native< MatrixS<1,3,MatrixS33>, VectorS<3,Real>, MatrixS33 >;
template class MatMul_Native< MatrixS<1,2,MatrixS44>, VectorS<2,Real>, MatrixS44 >;
template class MatMul_Native< MatrixS<1,2,MatrixS55>, VectorS<2,Real>, MatrixS55 >;
template class MatMul_Native< MatrixS<1,3,MatrixS55>, VectorS<3,Real>, MatrixS55 >;
template class MatMul_Native< MatrixS<1,3,MatrixS66>, VectorS<3,Real>, MatrixS66 >;
template class MatMul_Native< MatrixS<1,3,MatrixS77>, VectorS<3,Real>, MatrixS77 >;
// --

template class MatMul_Native< MatrixS<1,1,MatrixS22>, VectorS<1,MatrixS22>, MatrixS22 >;
template class MatMul_Native< MatrixS<1,2,MatrixS22>, VectorS<2,MatrixS22>, MatrixS22 >;
template class MatMul_Native< MatrixS<1,2,MatrixS44>, VectorS<2,MatrixS44>, MatrixS44 >;
template class MatMul_Native< MatrixS<1,2,MatrixS55>, VectorS<2,MatrixS55>, MatrixS55 >;
template class MatMul_Native< MatrixS<1,2,MatrixS66>, VectorS<2,MatrixS66>, MatrixS66 >;

template class MatMul_Native< MatrixS<1,3,MatrixS22>, VectorS<3,MatrixS22>, MatrixS22 >;
template class MatMul_Native< MatrixS<1,3,MatrixS33>, VectorS<3,MatrixS33>, MatrixS33 >;
template class MatMul_Native< MatrixS<1,3,MatrixS44>, VectorS<3,MatrixS44>, MatrixS44 >;
template class MatMul_Native< MatrixS<1,3,MatrixS55>, VectorS<3,MatrixS55>, MatrixS55 >;
template class MatMul_Native< MatrixS<1,3,MatrixS66>, VectorS<3,MatrixS66>, MatrixS66 >;
template class MatMul_Native< MatrixS<1,3,MatrixS77>, VectorS<3,MatrixS77>, MatrixS77 >;

template class MatMul_Native< MatrixS<2,2,MatrixS22>, VectorS<2,MatrixS22>, VectorS<2,MatrixS22> >;
template class MatMul_Native< MatrixS<2,2,MatrixS33>, VectorS<2,MatrixS33>, VectorS<2,MatrixS33> >;
template class MatMul_Native< MatrixS<2,2,MatrixS44>, VectorS<2,MatrixS44>, VectorS<2,MatrixS44> >;
template class MatMul_Native< MatrixS<2,2,MatrixS55>, VectorS<2,MatrixS55>, VectorS<2,MatrixS55> >;

template class MatMul_Native< MatrixS<3,3,MatrixS22>, VectorS<3,MatrixS22>, VectorS<3,MatrixS22> >;
template class MatMul_Native< MatrixS<3,3,MatrixS33>, VectorS<3,MatrixS33>, VectorS<3,MatrixS33> >;
template class MatMul_Native< MatrixS<3,3,MatrixS44>, VectorS<3,MatrixS44>, VectorS<3,MatrixS44> >;

template class MatMul_Native< MatrixS<1,1,MatrixS44>, VectorS<1,VectorS4>, VectorS4 >;
template class MatMul_Native< MatrixS<1,2,MatrixS44>, VectorS<2,VectorS4>, VectorS4 >;
template class MatMul_Native< MatrixS<1,2,MatrixS55>, VectorS<2,VectorS5>, VectorS5 >;
template class MatMul_Native< MatrixS<1,2,MatrixS66>, VectorS<2,VectorS6>, VectorS6 >;

template class MatMul_Native< MatrixS<1,2,MatrixS32>, VectorS<2,VectorS2>, VectorS<1,VectorS3> >;

template class MatMul_Native< MatrixS<1,3,MatrixS55>, VectorS<3,VectorS5>, VectorS5 >;
template class MatMul_Native< MatrixS<1,3,MatrixS66>, VectorS<3,VectorS6>, VectorS6 >;
template class MatMul_Native< MatrixS<1,3,MatrixS77>, VectorS<3,VectorS7>, VectorS7 >;

template class MatMul_Native< MatrixS<2,2,MatrixS22>, VectorS<2,VectorS2>, VectorS<2,VectorS2> >;
template class MatMul_Native< MatrixS<2,2,MatrixS33>, VectorS<2,VectorS3>, VectorS<2,VectorS3> >;
template class MatMul_Native< MatrixS<2,2,MatrixS44>, VectorS<2,VectorS4>, VectorS<2,VectorS4> >;
template class MatMul_Native< MatrixS<2,2,MatrixS55>, VectorS<2,VectorS5>, VectorS<2,VectorS5> >;


template class MatMul_Native< MatrixS<3,3,MatrixS22>, VectorS<3,VectorS2>, VectorS<3,VectorS2> >;
template class MatMul_Native< MatrixS<3,3,MatrixS33>, VectorS<3,VectorS3>, VectorS<3,VectorS3> >;
template class MatMul_Native< MatrixS<3,3,MatrixS44>, VectorS<3,VectorS4>, VectorS<3,VectorS4> >;

template class MatMul_Native< VectorS<1,Real>, MatrixS<1,1,Real>, VectorS<1,Real> >;
template class MatMul_Native< VectorS<2,Real>, MatrixS<1,1,Real>, VectorS<2,Real> >;
template class MatMul_Native< VectorS<3,Real>, MatrixS<1,1,Real>, VectorS<3,Real> >;
template class MatMul_Native< VectorS<4,Real>, MatrixS<1,1,Real>, VectorS<4,Real> >;

template class MatMul_Native< VectorS<1,MatrixS22>, MatrixS<1,1,VectorS2>, VectorS<1,VectorS2> >;
template class MatMul_Native< VectorS<1,MatrixS44>, MatrixS<1,1,VectorS4>, VectorS<1,VectorS4> >;
template class MatMul_Native< VectorS<2,MatrixS22>, MatrixS<1,1,VectorS2>, VectorS<2,VectorS2> >;
template class MatMul_Native< VectorS<2,MatrixS33>, MatrixS<1,1,VectorS3>, VectorS<2,VectorS3> >;
template class MatMul_Native< VectorS<2,MatrixS44>, MatrixS<1,1,VectorS4>, VectorS<2,VectorS4> >;
template class MatMul_Native< VectorS<2,MatrixS55>, MatrixS<1,1,VectorS5>, VectorS<2,VectorS5> >;
template class MatMul_Native< VectorS<2,MatrixS66>, MatrixS<1,1,VectorS6>, VectorS<2,VectorS6> >;

template class MatMul_Native< VectorS<3,MatrixS22>, MatrixS<1,1,VectorS2>, VectorS<3,VectorS2> >;
template class MatMul_Native< VectorS<3,MatrixS33>, MatrixS<1,1,VectorS3>, VectorS<3,VectorS3> >;
template class MatMul_Native< VectorS<3,MatrixS55>, MatrixS<1,1,VectorS5>, VectorS<3,VectorS5> >;
template class MatMul_Native< VectorS<3,MatrixS66>, MatrixS<1,1,VectorS6>, VectorS<3,VectorS6> >;
template class MatMul_Native< VectorS<3,MatrixS77>, MatrixS<1,1,VectorS7>, VectorS<3,VectorS7> >;

template class MatMul_Native< VectorS<1,MatrixS33>, MatrixS<1,1,VectorS3>, VectorS<1,VectorS3> >;

template class MatMul_Native< MatrixS<1,2,MatrixS66>, VectorS<2,Real>, MatrixS66 >;

} //namespace DLA
} //namespace numpack 
