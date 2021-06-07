// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(MATRIXS_DECOMPOSE_LUP_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "MatrixS_Decompose_LUP.h"
#include "MatrixS.h"
#include "VectorS.h"
#include "tinymat/dense/InverseLUP.h"

//Perform a LU decomposition without pivoting and unit diagonal on U

namespace tinymat 
{
namespace DLA
{

namespace Fixed
{
//A pivoting class to distinguish between a matrix of matrices and other datatypes
template< class T >
struct pivoter
{
  template<int M>
  static int pivot(int j, MatrixS<M,M,T>& Matrix)
  {
     int ipivot = Matrix.max_row_in_col(j, j);

     Matrix.swap_rows(ipivot, j);

     return ipivot;
  }
};
#if 0
template<int M, int N, typename T >
struct pivoter< MatrixS<M,N,T> >
{
  //Don't do any pivoting when working with a matrix of matrices for now
  static inline void pivot(int, MatrixDView<T>&, MatrixDView<T>& ) {}
};
#endif
}

//-----------------------------------------------------------------------------
template< int M, class T >
void MatrixSLUP<M, T>::Decompose( MatrixS< M, M, T >& Matrix, VectorS<M, int>& P )
{
  //Perform the LU decomposition
  for (int j = 0; j < M-1; ++j)
  {
    //Perform the pivot
    P[j] = Fixed::pivoter<T>::pivot(j, Matrix);

    T invdiag = InverseLUP::Inverse( Matrix(j ,j) );

    //Scale the rows for the Upper matrix
    Matrix.scale_row(j, invdiag, j+1);

    for (int i = j+1; i < M; ++i)
    {
      T factor = -Matrix(i, j);

      Matrix.axpy_rows(j, i, factor, j+1);
    }
  }

  P[M-1] = 0;
}

//-----------------------------------------------------------------------------
template< int M, class T >
MatrixSLUP<M, T>::MatrixSLUP(const DLA::Identity& I) : MatrixFac(I)
{
  for (int n = 0; n < M-1; n++)
    P[n] = n;
  P[M-1] = 0;
}

} //namespace DLA
} //namespace tinymat 
