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
#include "LinearAlgebra/SparseLinAlg/ScalarMatrix_CRS_impl.h"

#define SCALARVECTOR_INSTANTIATE
#include "LinearAlgebra/SparseLinAlg/ScalarVector_impl.h"

#include "LinearAlgebra/BlockLinAlg/MatrixBlock_3x3.h"
#include "LinearAlgebra/BlockLinAlg/VectorBlock_3.h"

#include "LinearAlgebra/DenseLinAlg/DynamicSize/MatrixD.h"

namespace SANS
{
namespace SLA
{

template class UMFPACK<
    BLA::MatrixBlock_3x3< DLA::MatrixD<SparseMatrix_CRS<Real                  >>,
                          DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,3,Real>>>,
                          DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<1,2,Real>>>,

                          DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,1,Real>>>,
                          DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,3,Real>>>,
                          DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<3,2,Real>>>,

                          DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,1,Real>>>,
                          DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,3,Real>>>,
                          DLA::MatrixD<SparseMatrix_CRS<DLA::MatrixS<2,2,Real>>> > >;

} //namespace SLA
} //namespace SANS
