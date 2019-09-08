// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define SCALARMATRIX_CRS_INSTANTIATE
#include "MKL_PARDISOSolver_defines.h"
#include "numpack/SparseLinAlg/ScalarMatrix_CRS_impl.h"

namespace numpack 
{
namespace SLA
{
template class ScalarMatrix_CRS<SANS_MKL_PARDISO_INT>;

template ScalarMatrix_CRS<SANS_MKL_PARDISO_INT>::ScalarMatrix_CRS(const SparseMatrix_CRS<Real>&);
template ScalarMatrix_CRS<SANS_MKL_PARDISO_INT>::ScalarMatrix_CRS(const DLA::MatrixDView<SparseMatrix_CRS<Real> >&);
}
}
