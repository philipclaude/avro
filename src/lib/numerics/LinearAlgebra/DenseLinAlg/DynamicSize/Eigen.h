// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXD_EIGEN_H
#define MATRIXD_EIGEN_H

#include "MatrixD_Type.h"

//Eigen value/vector interface for DenseLinAlg.

namespace SANS
{
namespace DLA
{
  template< class T >
  inline void
  EigenValues( MatrixDView<T>& A, MatrixDView<T>& wr, MatrixDView<T>& wi );

  template< class T >
  inline void
  EigenVectors( MatrixDView<T>& A, MatrixDView<T>& vl, MatrixDView<T>& vr );

  template< class T >
  inline void
  EigenSystem( MatrixDView<T>& A, MatrixDView<T>& wr, MatrixDView<T>& wi, MatrixDView<T>& vl, MatrixDView<T>& vr );

} //namespace SANS
} //namespace DLA

#include "lapack/Eigen.h"

#endif //MATRIXD_EIGEN_H
