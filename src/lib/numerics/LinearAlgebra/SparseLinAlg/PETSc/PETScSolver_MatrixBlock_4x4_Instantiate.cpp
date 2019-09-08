// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define PETSc_INSTANTIATE
//#include "PETScSolver_backsolve.h"
//#include "PETScSolver_factorize.h"
#include "PETScSolver_impl.h"

#define SCALARMATRIX_CRS_INSTANTIATE
#include "LinearAlgebra/SparseLinAlg/ScalarMatrix_CRS_impl.h"

#define SCALARVECTOR_INSTANTIATE
#include "LinearAlgebra/SparseLinAlg/ScalarVector_impl.h"

#include "LinearAlgebra/BlockLinAlg/MatrixBlock_4x4.h"
#include "LinearAlgebra/BlockLinAlg/VectorBlock_4.h"


namespace SANS
{
namespace SLA
{

//-----------------------------------------------------------------------------
template< class Matrix_type >
void PETScSolver<Matrix_type>::init()
{
  SANS_DEVELOPER_EXCEPTION("Not implemented");
}

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
  SANS_DEVELOPER_EXCEPTION("Not implemented"); return LinearSolveStatus();
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

template class PETScSolver<
                 BLA::MatrixBlock_4x4<
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<8,8,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<8,1,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<8,2,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<8,1,Real> > >,

                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<1,8,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<1,2,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,

                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<2,8,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<2,1,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<2,2,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<2,1,Real> > >,

                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<1,8,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<1,2,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >
                                      >
                          >;

template class PETScSolver<
                 BLA::MatrixBlock_4x4<
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<8,8,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<8,1,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<8,3,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<8,1,Real> > >,

                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<1,8,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<1,3,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,

                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<3,8,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<3,1,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<3,3,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<3,1,Real> > >,

                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<1,8,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<1,3,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >
                                      >
                          >;

template class PETScSolver<
                 BLA::MatrixBlock_4x4<
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<8,8,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<8,1,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<8,4,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<8,1,Real> > >,

                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<1,8,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<1,4,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,

                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<4,8,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<4,1,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<4,4,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<4,1,Real> > >,

                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<1,8,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<1,4,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >
                                      >
                          >;

template class PETScSolver<
                 BLA::MatrixBlock_4x4<
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<8,8,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<8,1,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<8,6,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<8,1,Real> > >,

                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<1,8,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<1,6,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,

                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<6,8,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<6,1,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<6,6,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<6,1,Real> > >,

                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<1,8,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<DLA::MatrixS<1,6,Real> > >,
                   DLA::MatrixD<SLA::SparseMatrix_CRS<Real> >
                                      >
                          >;

} //namespace SLA
} //namespace SANS
