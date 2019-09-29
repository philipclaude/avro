// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXS_TRACE_H
#define MATRIXS_TRACE_H

#include "MatrixS_Type.h"
#include "numpack/types/PromoteSurreal.h"

namespace numpack 
{
namespace DLA
{

  inline Real
  tr(const Real& A)
  {
    return A;
  }

  template<int N>
  inline SurrealS<N>
  tr(const SurrealS<N>& A)
  {
    return A;
  }

  //Computes the trace of MatrixS
  template<int M, class T>
  inline T
  tr(const MatrixS<M, M, T>& A)
  {
    T sum = 0;
    for (int i = 0; i < M; i++)
      sum += A(i,i);
    return sum;
  }

  //Computes the trace of MatrixSymS
  template<int M, class T>
  inline T
  tr(const MatrixSymS<M, T>& A)
  {
    T sum = 0;
    for (int i = 0; i < M; i++)
      sum += A(i,i);
    return sum;
  }

  //Computes the trace of tr(R*S) where R and S are MatrixS
  template<int M, int N, class TL, class TR>
  inline typename promote_Surreal<TL,TR>::type
  tr(const OpMulS< MatrixS<M,N,TL>, MatrixS<N,M,TR> >& tree )
  {
    typedef typename promote_Surreal<TL,TR>::type T;

    const MatrixS<M,N,TL>& ML = tree.left();
    const MatrixS<N,M,TR>& MR = tree.right();

    T sum = 0;
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        sum += ML(i,j)*MR(j,i);

    return sum;
  }

  //Computes the trace of tr(R*S) where R and S are MatrixSymS
  template<int M, class TL, class TR>
  inline typename promote_Surreal<TL,TR>::type
  tr(const OpMulS< MatrixSymS<M,TL>, MatrixSymS<M,TR> >& tree )
  {
    typedef typename promote_Surreal<TL,TR>::type T;

    const MatrixSymS<M,TL>& ML = tree.left();
    const MatrixSymS<M,TR>& MR = tree.right();

    T sum = 0;
    for (int i = 0; i < M; i++)
      for (int j = 0; j < M; j++)
        sum += ML(i,j)*MR(i,j);

    return sum;
  }

} //namespace DLA
} //namespace numpack 

#endif //MATRIXS_TRACE_H
