// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXD_INVERSELQR_INSTANTIATE
#include "MatrixD_InverseQR_impl.h"

#include "MatrixD_TupleMatrix.h"

namespace tinymat 
{
namespace DLA
{

template void ApplyQTrans<Real, MatrixDView< Real > >( MatrixDView< Real >& A, const VectorD<Real>& tau, MatrixDView< Real >& B );

// Explicitly instantiate all datatypes used here
template struct MatrixDQRSolver< Real, MatrixDView< Real > >;
//template struct MatrixDQRSolver< Real, MatrixDTuple< MatrixDView< Real > > >;

template struct MatrixDQRSolver< Real, MatrixDView< VectorS<1, Real> > >;

}
}
