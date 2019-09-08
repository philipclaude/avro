// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "LinearAlgebra/DenseLinAlg/StaticSize/MatrixSymS.h"
#include "types/SurrealS.h"

#define MATRIXD_MATMUL_NATIVE_INSTANTIATE
#include "MatMul_Native_impl.h"

namespace SANS
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
template class MatMul_Native< Real, Real, Real >;

template class MatMul_Native< Real, VectorS<1,Real>, VectorS<1,Real> >;
template class MatMul_Native< Real, VectorS<2,Real>, VectorS<2,Real> >;
template class MatMul_Native< Real, VectorS<3,Real>, VectorS<3,Real> >;
template class MatMul_Native< Real, VectorS<4,Real>, VectorS<4,Real> >;
template class MatMul_Native< Real, VectorS<5,Real>, VectorS<5,Real> >;
template class MatMul_Native< Real, VectorS<6,Real>, VectorS<6,Real> >;
template class MatMul_Native< Real, VectorS<7,Real>, VectorS<7,Real> >;
template class MatMul_Native< Real, VectorS<8,Real>, VectorS<8,Real> >;
template class MatMul_Native< Real, VectorS<22,Real>, VectorS<22,Real> >;

template class MatMul_Native< Real, VectorS<8,SurrealS<1>>, VectorS<8,SurrealS<1>> >;

template class MatMul_Native< Real, MatrixS<1,2,Real>, MatrixS<1,2,Real> >;
template class MatMul_Native< Real, MatrixS<1,3,Real>, MatrixS<1,3,Real> >;
template class MatMul_Native< Real, MatrixS<1,4,Real>, MatrixS<1,4,Real> >;
template class MatMul_Native< Real, MatrixS<1,5,Real>, MatrixS<1,5,Real> >;
template class MatMul_Native< Real, MatrixS<1,6,Real>, MatrixS<1,6,Real> >;
template class MatMul_Native< Real, MatrixS<1,7,Real>, MatrixS<1,7,Real> >;

template class MatMul_Native< Real, MatrixS<3,3,Real>, MatrixS<3,3,Real> >;

template class MatMul_Native< Real, MatrixSymS<1,Real>, MatrixSymS<1,Real> >;
template class MatMul_Native< Real, MatrixSymS<2,Real>, MatrixSymS<2,Real> >;
template class MatMul_Native< Real, MatrixSymS<3,Real>, MatrixSymS<3,Real> >;

template class MatMul_Native< Real, MatrixSymS<4,Real>, MatrixSymS<4,Real> >;

template class MatMul_Native< Real, SurrealS<1>, SurrealS<1> >;
template class MatMul_Native< Real, SurrealS<2>, SurrealS<2> >;

template class MatMul_Native< Real, MatrixS<3,3, SurrealS<1> >, MatrixS<3,3, SurrealS<1> > >;

typedef VectorS<2,Real> VectorS2;
typedef VectorS<3,Real> VectorS3;
typedef VectorS<4,Real> VectorS4;
typedef VectorS<5,Real> VectorS5;
typedef VectorS<6,Real> VectorS6;
typedef VectorS<7,Real> VectorS7;

template class MatMul_Native< Real, VectorS<1,VectorS2>, VectorS<1,VectorS2> >;
template class MatMul_Native< Real, VectorS<1,VectorS3>, VectorS<1,VectorS3> >;
template class MatMul_Native< Real, VectorS<1,VectorS4>, VectorS<1,VectorS4> >;

template class MatMul_Native< Real, VectorS<2,VectorS2>, VectorS<2,VectorS2> >;
template class MatMul_Native< Real, VectorS<2,VectorS3>, VectorS<2,VectorS3> >;
template class MatMul_Native< Real, VectorS<2,VectorS4>, VectorS<2,VectorS4> >;
template class MatMul_Native< Real, VectorS<2,VectorS5>, VectorS<2,VectorS5> >;
template class MatMul_Native< Real, VectorS<2,VectorS6>, VectorS<2,VectorS6> >;

template class MatMul_Native< Real, VectorS<3,VectorS2>, VectorS<3,VectorS2> >;
template class MatMul_Native< Real, VectorS<3,VectorS3>, VectorS<3,VectorS3> >;
template class MatMul_Native< Real, VectorS<3,VectorS4>, VectorS<3,VectorS4> >;
template class MatMul_Native< Real, VectorS<3,VectorS5>, VectorS<3,VectorS5> >;
template class MatMul_Native< Real, VectorS<3,VectorS6>, VectorS<3,VectorS6> >;
template class MatMul_Native< Real, VectorS<3,VectorS7>, VectorS<3,VectorS7> >;

template class MatMul_Native< Real, VectorS<4,VectorS2> , VectorS<4,VectorS2> >;

typedef MatrixS<2,2,Real> MatrixS22;
typedef MatrixS<3,3,Real> MatrixS33;
typedef MatrixS<4,4,Real> MatrixS44;
typedef MatrixS<5,5,Real> MatrixS55;
typedef MatrixS<6,6,Real> MatrixS66;

template class MatMul_Native< Real, VectorS<1,MatrixS33>, VectorS<1,MatrixS33> >;

template class MatMul_Native< Real, VectorS<2,MatrixS22>, VectorS<2,MatrixS22> >;
template class MatMul_Native< Real, VectorS<2,MatrixS33>, VectorS<2,MatrixS33> >;
template class MatMul_Native< Real, VectorS<2,MatrixS44>, VectorS<2,MatrixS44> >;
template class MatMul_Native< Real, VectorS<2,MatrixS55>, VectorS<2,MatrixS55> >;
template class MatMul_Native< Real, VectorS<2,MatrixS66>, VectorS<2,MatrixS66> >;

template class MatMul_Native< Real, VectorS<3,MatrixS22>, VectorS<3,MatrixS22> >;
template class MatMul_Native< Real, VectorS<3,MatrixS44>, VectorS<3,MatrixS44> >;

} //namespace DLA
} //namespace SANS
