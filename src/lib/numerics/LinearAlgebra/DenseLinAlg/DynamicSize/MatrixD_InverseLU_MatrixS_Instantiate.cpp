// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXD_INVERSELU_INSTANTIATE
#include "MatrixD_InverseLU_impl.h"

#include "MatrixD_TupleMatrix.h"

namespace SANS
{
namespace DLA
{

template struct MatrixDLUSolver< MatrixS<1,1,Real>, MatrixDView< VectorS<1,Real> > >;
template struct MatrixDLUSolver< MatrixS<2,2,Real>, MatrixDView< VectorS<2,Real> > >;
template struct MatrixDLUSolver< MatrixS<3,3,Real>, MatrixDView< VectorS<3,Real> > >;
template struct MatrixDLUSolver< MatrixS<4,4,Real>, MatrixDView< VectorS<4,Real> > >;
template struct MatrixDLUSolver< MatrixS<5,5,Real>, MatrixDView< VectorS<5,Real> > >;
template struct MatrixDLUSolver< MatrixS<6,6,Real>, MatrixDView< VectorS<6,Real> > >;
template struct MatrixDLUSolver< MatrixS<7,7,Real>, MatrixDView< VectorS<7,Real> > >;

template struct MatrixDLUSolver< MatrixS<1,1,Real>, MatrixDView< MatrixS<1,1,Real> > >;
template struct MatrixDLUSolver< MatrixS<2,2,Real>, MatrixDView< MatrixS<2,2,Real> > >;
template struct MatrixDLUSolver< MatrixS<3,3,Real>, MatrixDView< MatrixS<3,3,Real> > >;
template struct MatrixDLUSolver< MatrixS<4,4,Real>, MatrixDView< MatrixS<4,4,Real> > >;
template struct MatrixDLUSolver< MatrixS<5,5,Real>, MatrixDView< MatrixS<5,5,Real> > >;
template struct MatrixDLUSolver< MatrixS<6,6,Real>, MatrixDView< MatrixS<6,6,Real> > >;
template struct MatrixDLUSolver< MatrixS<7,7,Real>, MatrixDView< MatrixS<7,7,Real> > >;

template struct MatrixDLUSolver< MatrixS<2,2,MatrixS<2,2,Real>>, MatrixDView< MatrixS<2,2,MatrixS<2,2,Real>> > >;
template struct MatrixDLUSolver< MatrixS<2,2,MatrixS<3,3,Real>>, MatrixDView< MatrixS<2,2,MatrixS<3,3,Real>> > >;
template struct MatrixDLUSolver< MatrixS<2,2,MatrixS<4,4,Real>>, MatrixDView< MatrixS<2,2,MatrixS<4,4,Real>> > >;
template struct MatrixDLUSolver< MatrixS<2,2,MatrixS<5,5,Real>>, MatrixDView< MatrixS<2,2,MatrixS<5,5,Real>> > >;

template struct MatrixDLUSolver< MatrixS<3,3,MatrixS<2,2,Real>>, MatrixDView< MatrixS<3,3,MatrixS<2,2,Real>> > >;
template struct MatrixDLUSolver< MatrixS<3,3,MatrixS<3,3,Real>>, MatrixDView< MatrixS<3,3,MatrixS<3,3,Real>> > >;
template struct MatrixDLUSolver< MatrixS<3,3,MatrixS<4,4,Real>>, MatrixDView< MatrixS<3,3,MatrixS<4,4,Real>> > >;
template struct MatrixDLUSolver< MatrixS<3,3,MatrixS<5,5,Real>>, MatrixDView< MatrixS<3,3,MatrixS<5,5,Real>> > >;
template struct MatrixDLUSolver< MatrixS<3,3,MatrixS<6,6,Real>>, MatrixDView< MatrixS<3,3,MatrixS<6,6,Real>> > >;

typedef VectorS<2,Real> VectorS2;
typedef VectorS<3,Real> VectorS3;
typedef VectorS<4,Real> VectorS4;
typedef VectorS<5,Real> VectorS5;
typedef VectorS<6,Real> VectorS6;
typedef VectorS<7,Real> VectorS7;
template struct MatrixDLUSolver< MatrixS<1,1,Real>, MatrixDView< VectorS<1,VectorS2> > >;
template struct MatrixDLUSolver< MatrixS<1,1,Real>, MatrixDView< VectorS<1,VectorS3> > >;
template struct MatrixDLUSolver< MatrixS<1,1,Real>, MatrixDView< VectorS<1,VectorS4> > >;

template struct MatrixDLUSolver< MatrixS<2,2,Real>, MatrixDView< VectorS<2,VectorS2> > >;
template struct MatrixDLUSolver< MatrixS<2,2,Real>, MatrixDView< VectorS<2,VectorS3> > >;
template struct MatrixDLUSolver< MatrixS<2,2,Real>, MatrixDView< VectorS<2,VectorS4> > >;
template struct MatrixDLUSolver< MatrixS<2,2,Real>, MatrixDView< VectorS<2,VectorS5> > >;
template struct MatrixDLUSolver< MatrixS<2,2,Real>, MatrixDView< VectorS<2,VectorS6> > >;

template struct MatrixDLUSolver< MatrixS<3,3,Real>, MatrixDView< VectorS<3,VectorS2> > >;
template struct MatrixDLUSolver< MatrixS<3,3,Real>, MatrixDView< VectorS<3,VectorS3> > >;
template struct MatrixDLUSolver< MatrixS<3,3,Real>, MatrixDView< VectorS<3,VectorS5> > >;
template struct MatrixDLUSolver< MatrixS<3,3,Real>, MatrixDView< VectorS<3,VectorS6> > >;
template struct MatrixDLUSolver< MatrixS<3,3,Real>, MatrixDView< VectorS<3,VectorS7> > >;

typedef MatrixS<2,2,Real> MatrixS22;
typedef MatrixS<3,3,Real> MatrixS33;
typedef MatrixS<4,4,Real> MatrixS44;
typedef MatrixS<5,5,Real> MatrixS55;
typedef MatrixS<6,6,Real> MatrixS66;
typedef MatrixS<7,7,Real> MatrixS77;
template struct MatrixDLUSolver< MatrixS<1,1,Real>, MatrixDView< VectorS<1,MatrixS22> > >;
template struct MatrixDLUSolver< MatrixS<1,1,Real>, MatrixDView< VectorS<1,MatrixS33> > >;
template struct MatrixDLUSolver< MatrixS<1,1,Real>, MatrixDView< VectorS<1,MatrixS44> > >;

template struct MatrixDLUSolver< MatrixS<2,2,Real>, MatrixDView< VectorS<2,MatrixS22> > >;
template struct MatrixDLUSolver< MatrixS<2,2,Real>, MatrixDView< VectorS<2,MatrixS33> > >;
template struct MatrixDLUSolver< MatrixS<2,2,Real>, MatrixDView< VectorS<2,MatrixS44> > >;
template struct MatrixDLUSolver< MatrixS<2,2,Real>, MatrixDView< VectorS<2,MatrixS55> > >;
template struct MatrixDLUSolver< MatrixS<2,2,Real>, MatrixDView< VectorS<2,MatrixS66> > >;

template struct MatrixDLUSolver< MatrixS<3,3,Real>, MatrixDView< VectorS<3,MatrixS22> > >;
template struct MatrixDLUSolver< MatrixS<3,3,Real>, MatrixDView< VectorS<3,MatrixS33> > >;
template struct MatrixDLUSolver< MatrixS<3,3,Real>, MatrixDView< VectorS<3,MatrixS55> > >;
template struct MatrixDLUSolver< MatrixS<3,3,Real>, MatrixDView< VectorS<3,MatrixS66> > >;
template struct MatrixDLUSolver< MatrixS<3,3,Real>, MatrixDView< VectorS<3,MatrixS77> > >;

}
}
