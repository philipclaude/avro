// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(MATRIXD_BLAS_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "numpack/dense/dynamic/MatrixD.h"
#include "numpack/dense/tools/dense_BLAS.h"
#include "tools/SANSException.h"
#include "MatMul_BLAS.h"

#ifdef DLA_BLAS


namespace numpack 
{
namespace DLA
{

  template<class T>
  void MatMul_BLAS<T>::value(const MatrixDView<T>& ML, const MatrixDView<T>& MR,
                             const T sgn, MatrixDView<T>& res )
  {
    const T *A = &ML(0,0);
    const T *B = &MR(0,0);
    T *C = &res(0,0);

    const int m = ML.m();

    SANS_ASSERT(ML.n() == MR.m());

    const int k = ML.n();
    const int n = MR.n();

    const int Astride = ML.stride();
    const int Bstride = MR.stride();
    const int Cstride = res.stride();

    SANS_ASSERT_MSG( C != A && C != B, "The LHS MatrixD may not also be used in a matrix multiplication on the RHS.");

    if (n == 1)
      GEMV( CblasRowMajor, CblasNoTrans, m, k, sgn, A, Astride, B, Bstride, 0., C, Cstride );
    else
      GEMM( CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, sgn, A, Astride, B, Bstride, 0., C, Cstride );
  }

  template<class T>
  void MatMul_BLAS<T>::plus(const MatrixDView<T>& ML, const MatrixDView<T>& MR,
                            const T sgn, MatrixDView<T>& res )
  {
    const T *A = &ML(0,0);
    const T *B = &MR(0,0);
    T *C = &res(0,0);

    const int m = ML.m();

    SANS_ASSERT(ML.n() == MR.m());

    const int k = ML.n();
    const int n = MR.n();

    const int Astride = ML.stride();
    const int Bstride = MR.stride();
    const int Cstride = res.stride();

    SANS_ASSERT_MSG( C != A && C != B, "The LHS MatrixD may not also be used in a matrix multiplication on the RHS.");

    if (n == 1)
      GEMV( CblasRowMajor, CblasNoTrans, m, k, sgn, A, Astride, B, Bstride, 1., C, Cstride );
    else
      GEMM( CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, sgn, A, Astride, B, Bstride, 1., C, Cstride );
  }

}
}

#endif
