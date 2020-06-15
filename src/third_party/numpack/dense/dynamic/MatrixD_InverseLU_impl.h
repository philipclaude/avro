// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(MATRIXD_INVERSELU_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "numpack/dense/static/MatrixS.h"
#include "numpack/dense/static/VectorS.h"
#include "numpack/dense/static/MatrixSymS.h"

#include "MatrixD_InverseLU.h"
#include "MatrixD_Decompose_LU.h"
#include "numpack/dense/InverseLU.h"

//This computes computes a matrix inverse using LU decomposition without pivoting

namespace numpack
{
namespace DLA
{

//-----------------------------------------------------------------------------
template< class T, class MatrixType >
void MatrixDLUSolver<T, MatrixType>::Solve( const FactorType& Factorized, const Real sgn, MatrixType& res )
{
  const MatrixD<T>& AFac = Factorized.MatrixFac;

  SANS_ASSERT( AFac.m() == AFac.n() );
  SANS_ASSERT( AFac.n() == res.m() );

  const int m = AFac.m();
  const int n = res.n();

  int i, j;

  //Forward solve (L has non-one diagonals)
  for (j = 0; j < m; ++j)
  {
    for (i = 0; i < j; ++i)
    {
      T factor = -AFac(j, i);
      res.axpy_rows(i, j, factor, 0, n);
    }

    T invdiag = InverseLU::Inverse(AFac(j ,j));
    res.scale_row(j, invdiag, 0, n);
  }

  //Backward solve (U has one diagonals, no scale needed)
  for (j = m-2; j >= 0; --j)
    for (i = m-1; i > j; --i)
    {
      T factor = -AFac(j, i);
      res.axpy_rows(i, j, factor, 0, n);
    }

  if ( sgn != 1 )
    res *= sgn;
}

} //namespace DLA
} //namespace numpack
