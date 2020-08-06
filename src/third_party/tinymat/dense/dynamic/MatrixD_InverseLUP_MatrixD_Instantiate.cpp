// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXD_INVERSELUP_INSTANTIATE
#include "MatrixD_InverseLUP_impl.h"

#include "MatrixD_TupleMatrix.h"

namespace tinymat 
{
namespace DLA
{

template struct MatrixDLUPSolver< MatrixD<Real>, MatrixDView< VectorD<Real> > >;
template struct MatrixDLUPSolver< MatrixD<Real>, MatrixDView< MatrixD<Real> > >;

template struct MatrixDLUPSolver< MatrixD<MatrixS<2,2,Real> >, MatrixDView< VectorD<VectorS<2,Real>> > >;
template struct MatrixDLUPSolver< MatrixD<MatrixS<3,3,Real> >, MatrixDView< VectorD<VectorS<3,Real>> > >;
template struct MatrixDLUPSolver< MatrixD<MatrixS<4,4,Real> >, MatrixDView< VectorD<VectorS<4,Real>> > >;
template struct MatrixDLUPSolver< MatrixD<MatrixS<5,5,Real> >, MatrixDView< VectorD<VectorS<5,Real>> > >;
template struct MatrixDLUPSolver< MatrixD<MatrixS<6,6,Real> >, MatrixDView< VectorD<VectorS<6,Real>> > >;
template struct MatrixDLUPSolver< MatrixD<MatrixS<7,7,Real> >, MatrixDView< VectorD<VectorS<7,Real>> > >;

template struct MatrixDLUPSolver< MatrixD<MatrixS<3,3,Real> >, MatrixDView< MatrixD<MatrixS<3,1,Real>> > >;
template struct MatrixDLUPSolver< MatrixD<MatrixS<3,3,Real> >, MatrixDView< MatrixD<MatrixS<3,3,Real>> > >;
template struct MatrixDLUPSolver< MatrixD<MatrixS<5,5,Real> >, MatrixDView< MatrixD<MatrixS<5,5,Real>> > >;
}
}
