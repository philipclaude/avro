// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATMUL_BLAS_H
#define MATMUL_BLAS_H

#include "LinearAlgebra/DenseLinAlg/DynamicSize/MatrixD_Type.h"

namespace SANS
{
namespace DLA
{


template<class T>
class MatMul_BLAS
{

public:
//-----------------------------------------------------------------------------
  static void value(const MatrixDView<T>& ML, const MatrixDView<T>& MR,
                    const T sgn, MatrixDView<T>& res );

//-----------------------------------------------------------------------------
  static void plus(const MatrixDView<T>& ML, const MatrixDView<T>& MR,
                   const T sgn, MatrixDView<T>& res );

};

} //namespace DLA
} //namespace SANS

#endif //MATMUL_BLAS_H
