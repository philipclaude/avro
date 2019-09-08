// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(MATRIXS_INVERSELU_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "MatrixS.h"

//#include <boost/mpl/assert.hpp>

#include "MatrixS_InverseLU.h"
#include "MatrixS_Decompose_LU.h"
#include "numpack/DenseLinAlg/InverseLU.h"

//This computes computes a matrix inverse using LU decomposition without pivoting

namespace numpack 
{
namespace DLA
{

//-----------------------------------------------------------------------------
template< int M, int N, class T, class MatrixType >
void MatrixSLUSolver<M, N, T, MatrixType>::Solve( const FactorType& Factorized, const Real sgn, MatrixType& res )
{
  //BOOST_MPL_ASSERT_RELATION( M, ==, MatrixType::M );
  static const int NR = MatrixType::N;

  const MatrixS<M,N,T>& AFac = Factorized.MatrixFac;

  int i, j;

  //Forward solve (L has non-one diagonals)
  for (j = 0; j < M; ++j)
  {
    for (i = 0; i < j; ++i)
    {
      T factor = -AFac(j, i);
      res.axpy_rows(i, j, factor, 0, NR);
    }

    T invdiag = InverseLU::Inverse( AFac(j ,j) );
    res.scale_row(j, invdiag, 0, NR);
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
