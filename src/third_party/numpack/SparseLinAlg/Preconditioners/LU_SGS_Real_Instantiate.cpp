// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define LU_SGS_INSTANTIATE
#include "LU_SGS_impl.h"
#include "numpack/SparseLinAlg/SparseMatrix_CRS.h"

namespace numpack 
{
namespace SLA
{


template class LU_SGS< SparseMatrix_CRS<Real> >;


} //namespace SLA
} //namespace numpack 
