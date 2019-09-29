// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXD_INVERSELUP_INSTANTIATE
#include "MatrixD_InverseLUP_impl.h"

#include "MatrixD_TupleMatrix.h"

namespace numpack 
{
namespace DLA
{

template struct MatrixDLUPSolver< MatrixS<1,1,Real>, MatrixDView< VectorS<1,Real> > >;
template struct MatrixDLUPSolver< MatrixS<2,2,Real>, MatrixDView< VectorS<2,Real> > >;
template struct MatrixDLUPSolver< MatrixS<4,4,Real>, MatrixDView< VectorS<4,Real> > >;
template struct MatrixDLUPSolver< MatrixS<5,5,Real>, MatrixDView< VectorS<5,Real> > >;

template struct MatrixDLUPSolver< MatrixS<1,1,Real>, MatrixDView< MatrixS<1,1,Real> > >;
template struct MatrixDLUPSolver< MatrixS<2,2,Real>, MatrixDView< MatrixS<2,2,Real> > >;
template struct MatrixDLUPSolver< MatrixS<4,4,Real>, MatrixDView< MatrixS<4,4,Real> > >;
template struct MatrixDLUPSolver< MatrixS<5,5,Real>, MatrixDView< MatrixS<5,5,Real> > >;

}
}
