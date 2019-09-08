// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MKL_PARDISO_INSTANTIATE
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

template class MKL_PARDISO< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixS<1,1,Real> > > >;
template class MKL_PARDISO< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixS<2,2,Real> > > >;
template class MKL_PARDISO< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixS<3,3,Real> > > >;
template class MKL_PARDISO< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixS<4,4,Real> > > >;
template class MKL_PARDISO< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixS<5,5,Real> > > >;
template class MKL_PARDISO< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixS<6,6,Real> > > >;
template class MKL_PARDISO< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixS<7,7,Real> > > >;
template class MKL_PARDISO< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixS<8,8,Real> > > >;

template ScalarMatrix_CRS<SANS_MKL_PARDISO_INT>::ScalarMatrix_CRS(const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<1,1,Real> > >&);
template ScalarMatrix_CRS<SANS_MKL_PARDISO_INT>::ScalarMatrix_CRS(const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<2,2,Real> > >&);
template ScalarMatrix_CRS<SANS_MKL_PARDISO_INT>::ScalarMatrix_CRS(const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<3,3,Real> > >&);
template ScalarMatrix_CRS<SANS_MKL_PARDISO_INT>::ScalarMatrix_CRS(const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<4,4,Real> > >&);
template ScalarMatrix_CRS<SANS_MKL_PARDISO_INT>::ScalarMatrix_CRS(const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<5,5,Real> > >&);
template ScalarMatrix_CRS<SANS_MKL_PARDISO_INT>::ScalarMatrix_CRS(const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<6,6,Real> > >&);
template ScalarMatrix_CRS<SANS_MKL_PARDISO_INT>::ScalarMatrix_CRS(const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<7,7,Real> > >&);
template ScalarMatrix_CRS<SANS_MKL_PARDISO_INT>::ScalarMatrix_CRS(const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<8,8,Real> > >&);

} //namespace SLA
} //namespace SANS
