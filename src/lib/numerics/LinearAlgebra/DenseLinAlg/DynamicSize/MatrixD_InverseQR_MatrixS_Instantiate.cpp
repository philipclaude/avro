// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXD_INVERSELQR_INSTANTIATE
#include "MatrixD_InverseQR_impl.h"

//#include "MatrixD_TupleMatrix.h"

namespace SANS
{
namespace DLA
{

template struct MatrixDQRSolver< MatrixS<1,1,Real>, MatrixDView< VectorS<1,Real> > >;
template struct MatrixDQRSolver< MatrixS<2,2,Real>, MatrixDView< VectorS<2,Real> > >;
template struct MatrixDQRSolver< MatrixS<4,4,Real>, MatrixDView< VectorS<4,Real> > >;

template struct MatrixDQRSolver< MatrixS<1,1,Real>, MatrixDView< MatrixS<1,1,Real> > >;
template struct MatrixDQRSolver< MatrixS<2,2,Real>, MatrixDView< MatrixS<2,2,Real> > >;
template struct MatrixDQRSolver< MatrixS<4,4,Real>, MatrixDView< MatrixS<4,4,Real> > >;

}
}
