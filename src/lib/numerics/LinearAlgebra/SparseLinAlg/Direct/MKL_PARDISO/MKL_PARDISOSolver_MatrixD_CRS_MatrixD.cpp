// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MKL_PARDISO_INSTANTIATE
#include "MKL_PARDISOSolver_backsolve.h"
#include "MKL_PARDISOSolver_factorize.h"
#include "MKL_PARDISOSolver_impl.h"

namespace SANS
{
namespace SLA
{

// TODO: Non-zero pattern fix needed

//template class MKL_PARDISO< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixD<Real> > > >;

} //namespace SLA
} //namespace SANS
