// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXS_FROBNORM_H
#define MATRIXS_FROBNORM_H

#include "MatrixS_Type.h"

namespace numpack 
{
namespace DLA
{
  //Computes the Frobenius norm square of MatrixS
  template<int M, int N, class T>
  inline T
  FrobNormSq(const MatrixS<M, N, T>& A)
  {
    T sum = 0.0;
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        sum += A(i,j)*A(i,j);

    return sum;
  }

  //Computes the Frobenius norm of MatrixS
  template<int M, int N, class T>
  inline T
  FrobNorm(const MatrixS<M, N, T>& A) { return sqrt(FrobNormSq(A)); }

  //Computes the Frobenius norm of MatrixSymS
  template<int M, class T>
  inline T
  FrobNormSq(const MatrixSymS<M, T>& A)
  {
    T sum = 0.0;
    // for (int i = 0; i < M; i++)
    //   for (int j = 0; j < M; j++)
    //     sum += A(i,j)*A(i,j);

    // fewer flop version of the above, exploits the symmetry
    for (int i = 0; i < M; i++)
    {
      sum += A(i,i)*A(i,i); // diagonal
      for (int j = i+1; j < M; j++)
        sum += 2*A(i,j)*A(i,j); // above-diagonal counted twice
    }

    return sum;
  }

  //Computes the Frobenius norm of MatrixSymS
  template<int M, class T>
  inline T
  FrobNorm(const MatrixSymS<M, T>& A) { return sqrt(FrobNormSq(A)); }

} //namespace DLA
} //namespace numpack 

#endif //MATRIXS_FROBNORM_H
