// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(MATRIXS_DECOMPOSE_LU_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "MatrixS_Decompose_LU.h"
#include "MatrixS.h"
#include "numpack/DenseLinAlg/InverseLU.h"

//Perform a LU decomposition without pivoting and unit diagonal on U

namespace numpack 
{
namespace DLA
{

//-----------------------------------------------------------------------------
template< int M, class T >
void MatrixSLU<M, T>::Decompose( MatrixS< M, M, T >& Matrix )
{
  //Perform the LU decomposition
  for (int j = 0; j < M-1; ++j)
  {
    T invdiag = InverseLU::Inverse( Matrix(j ,j) );

    //Scale the rows for the Upper matrix
    Matrix.scale_row(j, invdiag, j+1);

    for (int i = j+1; i < M; ++i)
    {
      T factor = -Matrix(i, j);

      Matrix.axpy_rows(j, i, factor, j+1);
    }
  }
}

} //namespace DLA
} //namespace numpack 
