// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef LAPACK_EIGEN_H
#define LAPACK_EIGEN_H

#include "numpack/DenseLinAlg/DynamicSize/MatrixD_Type.h"

//Lapack eigen value/vector interface for DenseLinAlg.

namespace numpack 
{
namespace DLA
{
  template<class T>
  struct LAPACK_Eigen
  {
    static void Value( MatrixDView<T>& A, VectorDView<T>& wr, VectorDView<T>& wi );
    static void Vectors( MatrixDView<T>& A, MatrixDView<T>& vl, MatrixDView<T>& vr );
    static void System( MatrixDView<T>& A, VectorDView<T>& wr, VectorDView<T>& wi,
                                           MatrixDView<T>& vl, MatrixDView<T>& vr );
  };

  template< class T >
  inline void
  EigenValues( MatrixDView<T>& A, VectorDView<T>& wr, VectorDView<T>& wi )
  {
    LAPACK_Eigen<T>::Value(A, wr, wi);
  }

  template< class T >
  inline void
  EigenVectors( MatrixDView<T>& A, MatrixDView<T>& vl, MatrixDView<T>& vr )
  {
    LAPACK_Eigen<T>::Vectors(A, vl, vr);
  }

  template< class T >
  inline void
  EigenSystem( MatrixDView<T>& A, VectorDView<T>& wr, VectorDView<T>& wi, MatrixDView<T>& vl, MatrixDView<T>& vr )
  {
    LAPACK_Eigen<T>::System(A, wr, wi, vl, vr);
  }

} //namespace numpack 
} //namespace DLA

#endif //LAPACK_EIGEN_H
