// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MKL_PARDISO_INSTANTIATE
#include "MKL_PARDISOSolver_backsolve.h"
#include "MKL_PARDISOSolver_factorize.h"
#include "MKL_PARDISOSolver_impl.h"

#define SCALARMATRIX_CRS_INSTANTIATE
#include "numpack/sparse/ScalarMatrix_CRS_impl.h"

#define SCALARVECTOR_INSTANTIATE
#include "numpack/sparse/ScalarVector_impl.h"

#include "numpack/block/MatrixBlock_2x2.h"
#include "numpack/block/VectorBlock_2.h"


namespace numpack 
{
namespace SLA
{

template class MKL_PARDISO<
                BLA::MatrixBlock_2x2<
                  DLA::MatrixD<SparseMatrix_CRS<Real>>, DLA::MatrixD<SparseMatrix_CRS<Real>>,
                  DLA::MatrixD<SparseMatrix_CRS<Real>>, DLA::MatrixD<SparseMatrix_CRS<Real>> >  >;

template class MKL_PARDISO<
                BLA::MatrixBlock_2x2<
                  DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,2,Real>>>, DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,1,Real>>>,
                  DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,2,Real>>>, DLA::MatrixD<SparseMatrix_CRS<Real>>                   >  >;

template class MKL_PARDISO<
                BLA::MatrixBlock_2x2<
                  DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,2,Real>>>, DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,2,Real>>>,
                  DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,2,Real>>>, DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,2,Real>>> >  >;

template class MKL_PARDISO<
                 BLA::MatrixBlock_2x2<
                   DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<8,8,Real>>>, DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<8,8,Real>>>,
                   DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<8,8,Real>>>, DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<8,8,Real>>>
                                     >
                      >;


template class MKL_PARDISO<
                BLA::MatrixBlock_2x2<
                  DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,3,Real>>>, DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,2,Real>>>,
                  DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,3,Real>>>, DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,2,Real>>> >  >;

template class MKL_PARDISO<
                 BLA::MatrixBlock_2x2<
                   DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<8,8,Real>>>, DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<8,2,Real>>>,
                   DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,8,Real>>>, DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,2,Real>>>
                                     >
                      >;

template class MKL_PARDISO<
                BLA::MatrixBlock_2x2<
                  DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,3,Real>>>, DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,1,Real>>>,
                  DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,3,Real>>>, DLA::MatrixD<SparseMatrix_CRS<Real>> >  >;

template class MKL_PARDISO<
                BLA::MatrixBlock_2x2<
                  DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<4,4,Real>>>, DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<4,1,Real>>>,
                  DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,4,Real>>>, DLA::MatrixD<SparseMatrix_CRS<Real>> >  >;

} //namespace SLA
} //namespace numpack 
