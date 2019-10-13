// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define UMFPACK_INSTANTIATE
#include "UMFPACKSolver_backsolve.h"
#include "UMFPACKSolver_factorize.h"
#include "UMFPACKSolver_impl.h"
#include "UMFPACKSolver_MatrixD_Solve_impl.h"

#define SCALARMATRIX_CRS_INSTANTIATE
#include "numpack/sparse/ScalarMatrix_CRS_impl.h"

#define SCALARVECTOR_INSTANTIATE
#include "numpack/sparse/ScalarVector_impl.h"

namespace numpack 
{
namespace SLA
{
// TODO: Need to fix NonZeroPattern for SparseMatrix_CRS< DLA::MatrixD<Real> >

//template class UMFPACK< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixD< DLA::MatrixS<1,1,Real> > > > >;
//template class UMFPACK< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixD< DLA::MatrixS<2,2,Real> > > > >;
//template class UMFPACK< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixD< DLA::MatrixS<3,3,Real> > > > >;
//template class UMFPACK< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixD< DLA::MatrixS<4,4,Real> > > > >;

//template ScalarMatrix_CRS<SANS_UMFPACK_INT>::ScalarMatrix_CRS( const DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixD< DLA::MatrixS<1,1,Real> > > >&);
//template ScalarMatrix_CRS<SANS_UMFPACK_INT>::ScalarMatrix_CRS( const DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixD< DLA::MatrixS<2,2,Real> > > >&);
//template ScalarMatrix_CRS<SANS_UMFPACK_INT>::ScalarMatrix_CRS( const DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixD< DLA::MatrixS<3,3,Real> > > >&);
//template ScalarMatrix_CRS<SANS_UMFPACK_INT>::ScalarMatrix_CRS( const DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixD< DLA::MatrixS<4,4,Real> > > >&);

} //namespace SLA
} //namespace numpack 