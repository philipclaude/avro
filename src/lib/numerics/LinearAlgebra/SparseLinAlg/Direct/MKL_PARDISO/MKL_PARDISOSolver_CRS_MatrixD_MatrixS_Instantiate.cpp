// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MKL_PARDISO_INSTANTIATE
#include "LinearAlgebra/DenseLinAlg/StaticSize/MatrixS.h"
#include "LinearAlgebra/DenseLinAlg/DynamicSize/MatrixD.h"
#include "MKL_PARDISOSolver_backsolve.h"
#include "MKL_PARDISOSolver_factorize.h"
#include "MKL_PARDISOSolver_impl.h"

#define SCALARMATRIX_CRS_INSTANTIATE
#include "LinearAlgebra/SparseLinAlg/ScalarMatrix_CRS_impl.h"

#define SCALARVECTOR_INSTANTIATE
#include "LinearAlgebra/SparseLinAlg/ScalarVector_impl.h"

namespace SANS
{
namespace SLA
{
// TODO: Non-zero pattern fix needed

//template class MKL_PARDISO< SparseMatrix_CRS< DLA::MatrixD< DLA::MatrixS<1,1,Real> > > >;
//template class MKL_PARDISO< SparseMatrix_CRS< DLA::MatrixD< DLA::MatrixS<2,2,Real> > > >;
//template class MKL_PARDISO< SparseMatrix_CRS< DLA::MatrixD< DLA::MatrixS<4,4,Real> > > >;

//template ScalarMatrix_CRS<SANS_MKL_PARDISO_INT>::ScalarMatrix_CRS(const SparseMatrix_CRS< DLA::MatrixD< DLA::MatrixS<1,1,Real> > >&);
//template ScalarMatrix_CRS<SANS_MKL_PARDISO_INT>::ScalarMatrix_CRS(const SparseMatrix_CRS< DLA::MatrixD< DLA::MatrixS<2,2,Real> > >&);
//template ScalarMatrix_CRS<SANS_MKL_PARDISO_INT>::ScalarMatrix_CRS(const SparseMatrix_CRS< DLA::MatrixD< DLA::MatrixS<4,4,Real> > >&);

} //namespace SLA
} //namespace SANS
