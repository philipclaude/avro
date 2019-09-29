// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(MATRIXSYMS_DECOMPOSE_LDLT_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "MatrixSymS_Decompose_LDLT.h"
#include "numpack/dense/tools/SingularException.h"
#include "MatrixSymS.h"

//Perform a LDL^T decomposition on a symmetric matrix

// http://www.physics.arizona.edu/~restrepo/475A/Notes/sourcea-/node66.html
// https://en.wikipedia.org/wiki/Cholesky_decomposition#LDL_decomposition_2

namespace numpack 
{
namespace DLA
{

//-----------------------------------------------------------------------------
template< int M, class T >
void MatrixSymSLDLT<M, T>::Decompose( MatrixSymS<M, T>& A )
{
  VectorS<M,T> v;
  T s;

  //Perform the LDL^T decomposition
  for (int j = 0; j < M; j++)
  {
    for (int i = 0; i < j; i++)
      v[i] = A(j, i) * A(i, i);

    s = 0;
    for (int i = 0; i < j; i++)
      s += A(j, i) * v[i];

    A(j, j) = v[j] = A(j, j) - s;

    for (int k = j+1; k < M; k++)
    {
      s = 0;
      for (int i = 0; i < j; i++)
        s += A(k, i) * v[i];

      T numer = (A(k, j) - s);

      SANS_ASSERT_NONSINGULAR(v[j]);

      A(k, j) = numer/v[j];
    }
  }

}

} //namespace DLA
} //namespace numpack 
