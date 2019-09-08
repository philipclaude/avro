// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define UMFPACK_INSTANTIATE
#include "UMFPACKSolver_backsolve.h"
#include "UMFPACKSolver_factorize.h"
#include "UMFPACKSolver_impl.h"
#include "UMFPACKSolver_Solve_impl.h"

#define SCALARMATRIX_CRS_INSTANTIATE
#include "numpack/SparseLinAlg/ScalarMatrix_CRS_impl.h"

#define SCALARVECTOR_INSTANTIATE
#include "numpack/SparseLinAlg/ScalarVector_impl.h"

namespace numpack 
{
namespace SLA
{

template class UMFPACK< SparseMatrix_CRS< DLA::MatrixS<2,2,Real> > >;
template class UMFPACK< SparseMatrix_CRS< DLA::MatrixS<4,4,Real> > >;
template class UMFPACK< SparseMatrix_CRS< DLA::MatrixS<6,6,Real> > >;

template ScalarMatrix_CRS<SANS_UMFPACK_INT>::ScalarMatrix_CRS(const SparseMatrix_CRS< DLA::MatrixS<2,2,Real> >&);
template ScalarMatrix_CRS<SANS_UMFPACK_INT>::ScalarMatrix_CRS(const SparseMatrix_CRS< DLA::MatrixS<4,4,Real> >&);
template ScalarMatrix_CRS<SANS_UMFPACK_INT>::ScalarMatrix_CRS(const SparseMatrix_CRS< DLA::MatrixS<6,6,Real> >&);

} //namespace SLA
} //namespace numpack 
