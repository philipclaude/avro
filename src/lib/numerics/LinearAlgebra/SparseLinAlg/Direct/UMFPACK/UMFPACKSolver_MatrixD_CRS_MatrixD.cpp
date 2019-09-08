// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define UMFPACK_INSTANTIATE
#include "UMFPACKSolver_backsolve.h"
#include "UMFPACKSolver_factorize.h"
#include "UMFPACKSolver_impl.h"
#include "UMFPACKSolver_MatrixD_Solve_impl.h"

namespace SANS
{
namespace SLA
{

// TODO: Need to fix NonZeroPattern for SparseMatrix_CRS< DLA::MatrixD<Real> >

//template class UMFPACK< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixD<Real> > > >;

} //namespace SLA
} //namespace SANS
