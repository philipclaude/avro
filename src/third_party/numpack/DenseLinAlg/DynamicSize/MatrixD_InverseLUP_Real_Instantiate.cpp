// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXD_INVERSELUP_INSTANTIATE
#include "MatrixD_InverseLUP_impl.h"

namespace numpack 
{
namespace DLA
{

// Explicitly instantiate all datatypes used here
template struct MatrixDLUPSolver< Real, MatrixDView< Real > >;
template struct MatrixDLUPSolver< Real, MatrixDTuple< MatrixDView< Real > > >;
template struct MatrixDLUPSolver< Real, MatrixDTuple< MatrixDTuple< MatrixDView< VectorS<1, Real> > > > >;

template struct MatrixDLUPSolver< Real, MatrixDView< VectorS<1, Real> > >;
template struct MatrixDLUPSolver< Real, MatrixDView< VectorS<1, VectorS<1, Real> > > >;
template struct MatrixDLUPSolver< Real, MatrixDView< VectorS<2, Real> > >;
template struct MatrixDLUPSolver< Real, MatrixDView< VectorS<2, VectorS<1, Real> > > >;

template struct MatrixDLUPSolver< Real, MatrixDView< VectorS<3, Real> > >;

}
}
