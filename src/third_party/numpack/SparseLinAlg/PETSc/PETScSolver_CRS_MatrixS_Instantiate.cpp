// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define PETSc_INSTANTIATE
#include "PETScSolver_backsolve.h"
#include "PETScSolver_factorize.h"
#include "PETScSolver_impl.h"

#define SCALARMATRIX_CRS_INSTANTIATE
#include "numpack/SparseLinAlg/ScalarMatrix_CRS_impl.h"

#define SCALARVECTOR_INSTANTIATE
#include "numpack/SparseLinAlg/ScalarVector_impl.h"

namespace numpack 
{
namespace SLA
{

template class PETScSolver< SparseMatrix_CRS< DLA::MatrixS<2,2,Real> > >;
template class PETScSolver< SparseMatrix_CRS< DLA::MatrixS<4,4,Real> > >;
template class PETScSolver< SparseMatrix_CRS< DLA::MatrixS<5,5,Real> > >;
template class PETScSolver< SparseMatrix_CRS< DLA::MatrixS<6,6,Real> > >;

template ScalarMatrix_CRS<PetscInt>::ScalarMatrix_CRS(const SparseMatrix_CRS< DLA::MatrixS<2,2,Real> >&);
template ScalarMatrix_CRS<PetscInt>::ScalarMatrix_CRS(const SparseMatrix_CRS< DLA::MatrixS<4,4,Real> >&);
template ScalarMatrix_CRS<PetscInt>::ScalarMatrix_CRS(const SparseMatrix_CRS< DLA::MatrixS<5,5,Real> >&);
template ScalarMatrix_CRS<PetscInt>::ScalarMatrix_CRS(const SparseMatrix_CRS< DLA::MatrixS<6,6,Real> >&);

} //namespace SLA
} //namespace numpack 
