// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXS_SVD_LAPACK_INSTANTIATE
#define LAPACK_GESVD LAPACK_sgesvd

//Instantiate Lapack SVD routines for single precision

#ifdef DLA_LAPACK

#include "MatrixS_SVD_lapack_impl.h"

namespace tinymat 
{
namespace DLA
{

//1x1, 2x2, 3x3 are implemented locally, so instantiate other sizes
INSTANTIATE_SVD(2,3,float)
INSTANTIATE_SVD(3,2,float)

}
}

#endif