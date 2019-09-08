// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(MATRIXD_INVERSELUP_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "numpack/DenseLinAlg/StaticSize/MatrixS.h"
#include "numpack/DenseLinAlg/StaticSize/VectorS.h"

#include "MatrixD_InverseLUP.h"
#include "MatrixD_Decompose_LUP.h"
#include "MatrixD_TupleMatrix.h"
#include "numpack/DenseLinAlg/InverseLUP.h"

#include <memory>
#include <algorithm>

//This computes computes a matrix inverse using LU decomposition with pivoting

namespace numpack 
{
namespace DLA
{

//-----------------------------------------------------------------------------
template< class T, class MatrixType >
void MatrixDLUPSolver<T, MatrixType>::Solve(const FactorType& Factorized, const Real sgn, MatrixType& res )
{
  const MatrixD<T>& AFac = Factorized.MatrixFac;
  const VectorD<int>& P = Factorized.P;

  SANS_ASSERT( AFac.m() == AFac.n() );
  SANS_ASSERT( AFac.n() == res.m() );

  const int m = AFac.m();
  const int n = res.n();

  int i, j;
#if 0
  std::unique_ptr<int[]> maxpivot( new int[m] );

  if ( n > m )
    for (j = 0; j < m; ++j)
      maxpivot[j] = n-1;
  else
    for (j = 0; j < m; ++j)
      maxpivot[j] = std::min(j,n-1);
#endif

  //Perform the pivoting from LU decomposition
  for (j = 0; j < m-1; ++j)
  {
    //Perform the pivot
    int jpivot = P[j];
    res.swap_rows(jpivot, j);
//    if (jpivot > maxpivot[j])
//      for (k = 0; k < jpivot; ++k)
//        maxpivot[k] = std::min(jpivot,n-1);
  }

  //Forward solve (L has non-one diagonals)
  for (j = 0; j < m; ++j)
  {
    for (i = 0; i < j; ++i)
    {
      T factor = -AFac(j, i);
      res.axpy_rows(i, j, factor, 0, n); //maxpivot[i]+1); //+1 for end because it is related to size
    }

    T invdiag = InverseLUP::Inverse( AFac(j ,j) );
    res.scale_row(j, invdiag, 0, n); //maxpivot[j]+1); //+1 for end because it is related to size
  }

  //Backward solve (U has one diagonals, no scale needed)
  for (j = m-2; j >= 0; --j)
    for (i = m-1; i > j; --i)
    {
      T factor = -AFac(j, i);
      res.axpy_rows(i, j, factor);
    }

  if ( sgn != 1 )
    res *= sgn;
}

} //namespace DLA
} //namespace numpack 
