// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(MATRIXD_DECOMPOSE_LU_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "MatrixD_Decompose_LU.h"
#include "numpack/dense/dynamic/MatrixD.h"
#include "numpack/dense/static/MatrixS.h"
#include "numpack/dense/InverseLU.h"

//Perform a LU decomposition without pivoting and unit diagonal on U

namespace numpack 
{
namespace DLA
{

//-----------------------------------------------------------------------------
template< class T >
void MatrixDLU<T>::Decompose( MatrixDView< T >& Matrix )
{
  SANS_ASSERT( Matrix.m() == Matrix.n() );

  const int m = Matrix.m();

  //Perform the LU decomposition
  for (int j = 0; j < m-1; ++j)
  {
    T invdiag = InverseLU::Inverse(Matrix(j ,j));

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
