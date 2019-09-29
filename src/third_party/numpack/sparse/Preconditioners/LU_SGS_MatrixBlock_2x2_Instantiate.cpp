// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define LU_SGS_INSTANTIATE
#include "LU_SGS_impl.h"
#include "numpack/sparse/SparseMatrix_CRS.h"
#include "numpack/dense/static/VectorS.h"

namespace numpack 
{
namespace SLA
{

typedef BLA::MatrixBlock_2x2<DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >, DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,
                             DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >, DLA::MatrixD<SLA::SparseMatrix_CRS<Real> > > BlockMatrixRealReal;

template class LU_SGS<BlockMatrixRealReal>;

typedef DLA::MatrixS<8,8,Real> MatrixQ88;
typedef BLA::MatrixBlock_2x2<DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixQ88> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixQ88> >,
                             DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixQ88> >, DLA::MatrixD<SLA::SparseMatrix_CRS<MatrixQ88> > > BlockMatrix88;

template class LU_SGS<BlockMatrix88>;

} //namespace SLA
} //namespace numpack 
