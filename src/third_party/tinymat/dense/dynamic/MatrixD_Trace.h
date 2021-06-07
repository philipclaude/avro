// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXS_TRACE_H
#define MATRIXS_TRACE_H

#include "MatrixD_Type.h"
#include "tinymat/types/PromoteSurreal.h"

namespace tinymat
{
namespace DLA
{

  template<int N>
  inline SurrealS<N>
  tr(const SurrealS<N>& A)
  {
    return A;
  }

  //Computes the trace of MatrixS
  template<class T>
  inline T
  tr(const MatrixD<T>& A)
  {
    const int M = A.m();
    T sum = 0;
    for (int i = 0; i < M; i++)
      sum += A(i,i);
    return sum;
  }

  //Computes the trace of MatrixSymS
  template<class T>
  inline T
  tr(const MatrixSymD<T>& A)
  {
    const int M = A.m();
    T sum = 0;
    for (int i = 0; i < M; i++)
      sum += A(i,i);
    return sum;
  }

  //Computes the trace of tr(R*S) where R and S are MatrixS
  template<class TL, class TR>
  inline typename promote_Surreal<TL,TR>::type
  tr(const OpMulS< MatrixD<TL>, MatrixD<TR> >& tree )
  {
    typedef typename promote_Surreal<TL,TR>::type T;

    //const MatrixS<M,N,TL>& ML = tree.left();
    //const MatrixS<N,M,TR>& MR = tree.right();

    const MatrixD<TL>& ML = tree.left();
    const MatrixD<TR>& MR = tree.right();

    const int M = ML.m();
    const int N = ML.n();

    assert( M == MR.n() );
    assert( N == MR.m() );

    T sum = 0;
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        sum += ML(i,j)*MR(j,i);

    return sum;
  }

  //Computes the trace of tr(R*S) where R and S are MatrixSymS
  template<class TL, class TR>
  inline typename promote_Surreal<TL,TR>::type
  tr(const OpMulS< MatrixSymD<TL>, MatrixSymD<TR> >& tree )
  {
    typedef typename promote_Surreal<TL,TR>::type T;

    const MatrixSymD<TL>& ML = tree.left();
    const MatrixSymD<TR>& MR = tree.right();

    const int M = ML.m();
    assert( MR.m() == M );

    T sum = 0;
    for (int i = 0; i < M; i++)
      for (int j = 0; j < M; j++)
        sum += ML(i,j)*MR(i,j);

    return sum;
  }

} //namespace DLA
} //namespace tinymat

#endif //MATRIXS_TRACE_H
