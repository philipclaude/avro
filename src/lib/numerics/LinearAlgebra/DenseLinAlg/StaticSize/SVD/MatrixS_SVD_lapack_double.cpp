// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXS_SVD_LAPACK_INSTANTIATE
#define LAPACK_GESVD LAPACK_dgesvd

//Instantiate Lapack SVD routines for double precision

#ifdef DLA_LAPACK

#include "MatrixS_SVD_lapack_impl.h"

namespace SANS
{
namespace DLA
{

//1x1, 2x2, 3x3 are implemented locally, so instantiate other sizes
INSTANTIATE_SVD(2,3,Real)
INSTANTIATE_SVD(3,2,Real)
INSTANTIATE_SVD(3,4,Real)
INSTANTIATE_SVD(4,4,Real)

}
}

#endif
