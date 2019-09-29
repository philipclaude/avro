// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "Matrix_Util.h"

#include "dense_BLAS.h"

#ifdef DLA_BLAS

namespace numpack 
{
namespace DLA
{

//=============================================================================
// BLAS float implementation
void MatrixUtil_BLAS<float>::swap(float* __restrict x, float* __restrict y, const int size)
{
  cblas_sswap(size, x, 1, y, 1);
}

void MatrixUtil_BLAS<float>::scal(float* __restrict x, const float& a, const int size)
{
  cblas_sscal(size, a, x, 1);
}

void MatrixUtil_BLAS<float>::axpy(const float* __restrict x, float* __restrict y, const float& a, const int size)
{
  cblas_saxpy(size, a, x, 1, y, 1);
}

int MatrixUtil_BLAS<float>::max_row_in_col(const float* __restrict x, const int size, const int stride, const int start)
{
  return cblas_isamax(size-start, x + start*stride, stride) + start;
}


//=============================================================================
// BLAS double implementation
void MatrixUtil_BLAS<double>::swap(double* __restrict x, double* __restrict y, const int size)
{
  cblas_dswap(size, x, 1, y, 1);
}

void MatrixUtil_BLAS<double>::scal(double* __restrict x, const double& a, const int size)
{
  cblas_dscal(size, a, x, 1);
}

void MatrixUtil_BLAS<double>::axpy(const double* __restrict x, double* __restrict y, const double& a, const int size)
{
  cblas_daxpy(size, a, x, 1, y, 1);
}

int MatrixUtil_BLAS<double>::max_row_in_col(const double* __restrict x, const int size, const int stride, const int start)
{
  return cblas_idamax(size-start, x + start*stride, stride) + start;
}

} //namespace DLA
} //namespace numpack 

#endif //DLA_BLAS
