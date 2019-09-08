// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXD_BLAS_INSTANTIATE
#define GEMV cblas_dgemv
#define GEMM cblas_dgemm

#include "LinearAlgebra/DenseLinAlg/tools/DenseLinAlg_BLAS.h"

#ifdef DLA_BLAS

#include "MatMul_BLAS_impl.h"

namespace SANS
{
namespace DLA
{
template class MatMul_BLAS<double>;
}
}

#endif
