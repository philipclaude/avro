// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXD_INVERSELU_INSTANTIATE
#include "MatrixD_InverseLU_impl.h"

#include "VectorD.h"
#include "MatrixD_TupleMatrix.h"

namespace SANS
{
namespace DLA
{

template struct MatrixDLUSolver< MatrixD<Real>, MatrixDView< VectorD<Real> > >;
template struct MatrixDLUSolver< MatrixD<Real>, MatrixDView< MatrixD<Real> > >;

template struct MatrixDLUSolver< Real, MatrixDView< VectorS<1, VectorS<4,Real>> > >;
template struct MatrixDLUSolver< Real, MatrixDView< VectorS<2, VectorS<3,Real>> > >;
template struct MatrixDLUSolver< Real, MatrixDView< VectorS<2, VectorS<4,Real>> > >;
template struct MatrixDLUSolver< Real, MatrixDView< VectorS<3, VectorS<2,Real>> > >;
template struct MatrixDLUSolver< Real, MatrixDView< VectorS<3, VectorS<4,Real>> > >;

template struct MatrixDLUSolver< MatrixD<MatrixS<2,2,Real> >, MatrixDView< VectorD<VectorS<2,Real>> > >;
template struct MatrixDLUSolver< MatrixD<MatrixS<3,3,Real> >, MatrixDView< VectorD<VectorS<3,Real>> > >;
template struct MatrixDLUSolver< MatrixD<MatrixS<4,4,Real> >, MatrixDView< VectorD<VectorS<4,Real>> > >;
template struct MatrixDLUSolver< MatrixD<MatrixS<5,5,Real> >, MatrixDView< VectorD<VectorS<5,Real>> > >;
template struct MatrixDLUSolver< MatrixD<MatrixS<6,6,Real> >, MatrixDView< VectorD<VectorS<6,Real>> > >;

}
}
