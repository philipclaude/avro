// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "MatrixD_Mul.h"
#include "MatrixD.h"

#ifdef DLA_BLAS

namespace tinymat 
{
namespace DLA
{

//-----------------------------------------------------------------------------
template<>
void OpMulD_impl<Real, Real, Real>::value(const MatrixDView<Real>& ML, const MatrixDView<Real>& MR, const Real sgn, MatrixDView<Real>& res )
{
  const int mL = ML.m();
  const int nL = ML.n();
  const int nR = MR.n();

  if ( nL == 0 )
    res = 0;
  else if (mL*nR > 1) //Could add a lower limit for the size to call BLAS here
    MatMul_BLAS<Real>::value(ML, MR, sgn, res);
  else
    MatMul_Native<Real, Real, Real>::value(ML, MR, sgn, res);
}

//-----------------------------------------------------------------------------
template<>
void OpMulD_impl<Real, Real, Real>::plus(const MatrixDView<Real>& ML, const MatrixDView<Real>& MR, const Real sgn, MatrixDView<Real>& res )
{
  const int mL = ML.m();
  const int nL = ML.n();
  const int nR = MR.n();

  if ( nL == 0 )
    return;
  else if (mL*nR > 1) //Could add a lower limit for the size to call BLAS here
    MatMul_BLAS<Real>::plus(ML, MR, sgn, res);
  else
    MatMul_Native<Real, Real, Real>::plus(ML, MR, sgn, res);
}

} //namespace DLA
} //namespace tinymat 
#endif
