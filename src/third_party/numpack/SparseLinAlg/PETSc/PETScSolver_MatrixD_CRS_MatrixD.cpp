// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "numpack/SparseLinAlg/SparseMatrix_CRS.h"
#include "numpack/SparseLinAlg/ScalarMatrix_CRS.h"


#define PETSc_INSTANTIATE
//#include "PETScSolver_backsolve.h"
//#include "PETScSolver_factorize.h"
#include "PETScSolver_impl.h"


namespace numpack 
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
void PETScSolver<Matrix_type>::factorizeMatrix()
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

template class PETScSolver< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixD<Real> > > >;

} //namespace SLA
} //namespace numpack 
