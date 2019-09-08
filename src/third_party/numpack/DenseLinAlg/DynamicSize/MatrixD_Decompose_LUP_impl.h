// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(MATRIXD_DECOMPOSE_LUP_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "MatrixD_Decompose_LUP.h"
#include "numpack/DenseLinAlg/StaticSize/MatrixS.h"
#include "numpack/DenseLinAlg/InverseLUP.h"

//LU decomposition with pivoting and a unit diagonal on U

namespace numpack 
{
namespace DLA
{

//A pivoting class to distinguish between a matrix of matrices and other datatypes
template< class T >
struct pivoter
{
  static int pivot(int j, MatrixDView<T>& Matrix)
  {
     int ipivot = Matrix.max_row_in_col(j, j);

     Matrix.swap_rows(ipivot, j);

     return ipivot;
  }
};

template<int M, int N, typename T >
struct pivoter< MatrixS<M,N,T> >
{
  //Don't do any pivoting when working with a matrix of matrices for now
  static inline int pivot(int j, MatrixDView< MatrixS<M,N,T> >&) { return j; }
};

template<typename T >
struct pivoter< MatrixD<T> >
{
  //Don't do any pivoting when working with a matrix of matrices for now
  static inline int pivot(int j, MatrixDView< MatrixD<T> >&) { return j; }
};

//-----------------------------------------------------------------------------
template< class T >
void MatrixDLUP<T>::Decompose( MatrixDView< T >& Matrix, VectorDView<int>& P )
{
  SANS_ASSERT( Matrix.m() == Matrix.n() );
  SANS_ASSERT( Matrix.m() == P.m() );

  const int m = Matrix.m();

  for (int i = 0; i < m; ++i)
    P[i] = i;

  //Perform the LU decomposition
  for (int j = 0; j < m-1; ++j)
  {
    //Perform the pivot
    int ipivot = pivoter<T>::pivot(j, Matrix);
    P[ipivot] = P[j];
    P[j] = ipivot;

    T invdiag = InverseLUP::Inverse( Matrix(j ,j) );

    //Scale the rows for the Upper matrix
    Matrix.scale_row(j, invdiag, j+1);

    for (int i = j+1; i < m; ++i)
    {
      T factor = -Matrix(i, j);

      Matrix.axpy_rows(j, i, factor, j+1);
    }
  }
}

} //namespace DLA
} //namespace numpack 
