// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXD_DECOMPOSE_QR_INSTANTIATE
#include "MatrixD_Decompose_QR_impl.h"

namespace numpack 
{
namespace DLA
{

template struct MatrixDQR<Real>;

}
}