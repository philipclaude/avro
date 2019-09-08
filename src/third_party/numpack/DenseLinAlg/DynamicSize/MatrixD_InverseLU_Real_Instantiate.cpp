// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXD_INVERSELU_INSTANTIATE
#include "MatrixD_InverseLU_impl.h"

#include "numpack/types/SurrealS.h"

namespace numpack 
{
namespace DLA
{

template struct MatrixDLUSolver< Real, MatrixDView< Real > >;
template struct MatrixDLUSolver< Real, MatrixDView< SurrealS<1> > >;
template struct MatrixDLUSolver< Real, MatrixDTuple< MatrixDView< Real > > >;
template struct MatrixDLUSolver< Real, MatrixDTuple< MatrixDTuple< MatrixDView< VectorS<1, Real> > > > >;

template struct MatrixDLUSolver< Real, MatrixDView< VectorS< 1, Real> > >;
template struct MatrixDLUSolver< Real, MatrixDView< VectorS< 2, Real> > >;
template struct MatrixDLUSolver< Real, MatrixDView< VectorS< 3, Real> > >;
template struct MatrixDLUSolver< Real, MatrixDView< VectorS< 4, Real> > >;
template struct MatrixDLUSolver< Real, MatrixDView< VectorS< 5, Real> > >;
template struct MatrixDLUSolver< Real, MatrixDView< VectorS< 6, Real> > >;
template struct MatrixDLUSolver< Real, MatrixDView< VectorS< 7, Real> > >;
template struct MatrixDLUSolver< Real, MatrixDView< VectorS< 8, Real> > >;
template struct MatrixDLUSolver< Real, MatrixDView< VectorS<10, Real> > >;
template struct MatrixDLUSolver< Real, MatrixDView< VectorS<22, Real> > >;

template struct MatrixDLUSolver< Real, MatrixDView< VectorS<1, VectorS<1, Real> > > >;
template struct MatrixDLUSolver< Real, MatrixDView< VectorS<1, VectorS<2, Real> > > >;
template struct MatrixDLUSolver< Real, MatrixDView< VectorS<1, VectorS<3, Real> > > >;
template struct MatrixDLUSolver< Real, MatrixDView< VectorS<2, VectorS<2, Real> > > >;
template struct MatrixDLUSolver< Real, MatrixDView< VectorS<2, VectorS<5, Real> > > >;
template struct MatrixDLUSolver< Real, MatrixDView< VectorS<2, VectorS<6, Real> > > >;
template struct MatrixDLUSolver< Real, MatrixDView< VectorS<3, VectorS<3, Real> > > >;
template struct MatrixDLUSolver< Real, MatrixDView< VectorS<3, VectorS<5, Real> > > >;
template struct MatrixDLUSolver< Real, MatrixDView< VectorS<3, VectorS<6, Real> > > >;
template struct MatrixDLUSolver< Real, MatrixDView< VectorS<3, VectorS<7, Real> > > >;
template struct MatrixDLUSolver< Real, MatrixDView< VectorS<4, VectorS<2, Real> > > >;

template struct MatrixDLUSolver< Real, MatrixDView< MatrixS<3, 3, Real> > >;
template struct MatrixDLUSolver< Real, MatrixDView< MatrixS<3, 3, SurrealS<1> > > >;

template struct MatrixDLUSolver< Real, MatrixDView< MatrixSymS<1, Real> > >;
template struct MatrixDLUSolver< Real, MatrixDView< MatrixSymS<2, Real> > >;
template struct MatrixDLUSolver< Real, MatrixDView< MatrixSymS<3, Real> > >;
template struct MatrixDLUSolver< Real, MatrixDView< MatrixSymS<4, Real> > >;

}
}
