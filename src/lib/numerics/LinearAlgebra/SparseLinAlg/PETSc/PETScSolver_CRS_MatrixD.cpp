// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define PETSc_INSTANTIATE
#include "PETScSolver_impl.h"
//#include "PETScSolver_factorize.h"
//#include "PETScSolver_backsolve.h"

#define SCALARMATRIX_CRS_INSTANTIATE
#include "LinearAlgebra/SparseLinAlg/ScalarMatrix_CRS_impl.h"

#define SCALARVECTOR_INSTANTIATE
#include "LinearAlgebra/SparseLinAlg/ScalarVector_impl.h"

#include "LinearAlgebra/DenseLinAlg/StaticSize/MatrixS.h"
#include "LinearAlgebra/DenseLinAlg/DynamicSize/MatrixD.h"

namespace SANS
{
namespace SLA
{

template< class Matrix_type >
LinearSolveStatus PETScSolver<Matrix_type>::solve(SparseVectorView_type& b, SparseVectorView_type& x)
{
  SANS_DEVELOPER_EXCEPTION("Not implemented");
  return LinearSolveStatus(false);
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
void PETScSolver<Matrix_type>::factorize( SparseVectorView_type& bcondensed, bool transpose )
{
  SANS_DEVELOPER_EXCEPTION("Not implemented");
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
LinearSolveStatus PETScSolver<Matrix_type>::
backsolve(const bool transpose, const SparseVectorView_type& b, SparseVectorView_type& x ) const
{
  SANS_DEVELOPER_EXCEPTION("Not implemented");

  return LinearSolveStatus();
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
void PETScSolver<Matrix_type>::factorize()
{
  SANS_DEVELOPER_EXCEPTION("Not implemented");
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
PetscErrorCode
PETScSolver<Matrix_type>::MatGetOrdering_MDF(Mat mat, MatOrderingType type, IS *irow, IS *icol)
{
  SANS_DEVELOPER_EXCEPTION("Not implemented");
  PetscFunctionReturn(0);
}

template class PETScSolver< SparseMatrix_CRS< DLA::MatrixD<Real> > >;

} //namespace SLA
} //namespace SANS
