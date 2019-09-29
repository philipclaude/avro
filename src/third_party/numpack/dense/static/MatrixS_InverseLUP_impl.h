// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(MATRIXS_INVERSELUP_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "MatrixS.h"
#include "VectorS.h"

//#include <boost/mpl/assert.hpp>

#include "MatrixS_InverseLUP.h"
#include "MatrixS_Decompose_LUP.h"
#include "numpack/dense/InverseLUP.h"

//This computes computes a matrix inverse using LU decomposition with pivoting

namespace numpack 
{
namespace DLA
{

//-----------------------------------------------------------------------------
template< int M, int N, class T, class MatrixType >
void MatrixSLUPSolver<M, N, T, MatrixType>::Solve( const FactorType& Factorized, const Real sgn, MatrixType& res )
{
  //BOOST_MPL_ASSERT_RELATION( M, ==, MatrixType::M );
  static const int NR = MatrixType::N;

  const MatrixS<M,N,T>& AFac = Factorized.MatrixFac;
  const VectorS<M,int>& P = Factorized.P;

  int i, j;
#if 0
  int maxpivot[M];

  if ( N > M )
    for (j = 0; j < M; ++j)
      maxpivot[j] = NR-1;
  else
    for (j = 0; j < M; ++j)
      maxpivot[j] = std::min(j,NR-1);
#endif

  //Perform the pivoting from LU decomposition
  for (j = 0; j < M-1; ++j)
  {
    //Perform the pivot
    int jpivot = P[j];
    res.swap_rows(jpivot, j);
#if 0
    if (jpivot > maxpivot[j])
      for (k = j; k < jpivot; ++k)
        maxpivot[k] = std::min(jpivot,NR-1);
#endif
  }

  //Forward solve (L has non-one diagonals)
  for (j = 0; j < M; ++j)
  {
    for (i = 0; i < j; ++i)
    {
      T factor = -AFac(j, i);
      res.axpy_rows(i, j, factor, 0, NR); //maxpivot[i]+1); //+1 for end because it is related to size
    }

    T invdiag = InverseLUP::Inverse( AFac(j ,j) );
    res.scale_row(j, invdiag, 0, NR); //maxpivot[j]+1); //+1 for end because it is related to size
  }

  //Backward solve (U has one diagonals, no scale needed)
  for (j = M-2; j >= 0; --j)
    for (i = M-1; i > j; --i)
    {
      T factor = -AFac(j, i);
      res.axpy_rows(i, j, factor, 0, NR);
    }

  if ( sgn != 1 )
    res *= sgn;
}

} //namespace DLA
} //namespace numpack 
