// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(PETSc_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "tools/timer.h"

#include "numpack/sparse/ScalarMatrix_CRS.h"
#include "numpack/sparse/WriteMatrixMarketFile.h"

#include "PETScSolver.h"
#include "PETSc_VectorSize_impl.h"
#include "MDFOrdering.h"

#include <set>
#include <limits>
#include <petscmat.h>
#include <petsc/private/matorderimpl.h>

#define MPI_COMMUNICATOR_IN_CPP
#include "MPI/communicator.h"

namespace numpack 
{
namespace SLA
{


template< class Matrix_type >
LinearSolveStatus PETScSolver<Matrix_type>::solve(SparseVectorView_type& b, SparseVectorView_type& x)
{
  // factorize the matrix and then solve
  factorize();
  return backsolve( b, x );
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
void PETScSolver<Matrix_type>::factorize( SparseVectorView_type& bcondensed, bool transpose )
{
  SANS_DEVELOPER_EXCEPTION("Static Condensation not applicable for this Matrix Type");
}


//-----------------------------------------------------------------------------
template< class Matrix_type >
void PETScSolver<Matrix_type>::init()
{
  SystemNonZeroPattern nz(f_.matrixSize());
  f_.jacobian(nz);

  SANS_ASSERT( continuousmap_.size() == 1); //for now...
  SANS_ASSERT( staticCondensed_ == false); //no static condensation unless we're doing matrixD's

  // construct the matrix
  this->A_ = new Matrix_type(nz);

  // Flatten to a scalar matrix
  pMs_ = new ScalarMatrix_CRS<PetscInt, PetscScalar>(*this->A_);

  PetscInt m = pMs_->m();
  PetscInt n = pMs_->n();

  SANS_ASSERT( m <= n );  // Can only solve square matrices (but the matrix may contain off processor ghost columns)

  PetscScalar *a = pMs_->Rx();
  PetscInt *row_ptr = pMs_->Rp();
  PetscInt *col_ind = pMs_->Ri();

  static const int M = StateVectorSize<SparseVectorView_type>::M;

  // convert the local column indexing to a global continuous index
  int nnz = pMs_->nnz();
  for (int i = 0; i < nnz; i++)
  {
    if ( col_ind[i] < m )
      col_ind[i] += M*continuousmap_[0].nDOF_rank_offset;
    else
      col_ind[i] = col_ind[i] % M + M*continuousmap_[0].remoteGhostIndex[col_ind[i]/M - continuousmap_[0].nDOFpossessed];
  }

  // Create the PETSc matrix
#ifdef SANS_MPI
  // PETSc copies the memory with this call
  PETSc_STATUS( MatCreateMPIAIJWithArrays(*continuousmap_[0].comm, m, m, PETSC_DETERMINE, PETSC_DETERMINE, row_ptr, col_ind, a, &A_petsc_) );
  delete pMs_; pMs_ = nullptr;
#else
  // PETSc uses the provided memory with this call
  PETSc_STATUS( MatCreateSeqAIJWithArrays(*continuousmap_[0].comm, m, n, row_ptr, col_ind, a, &A_petsc_) );
#endif

  init_petsc();
}


//-----------------------------------------------------------------------------
template< class Matrix_type >
void PETScSolver<Matrix_type>::factorize()
{
  SANS_ASSERT( staticCondensed_ == false );
  factorizeMatrix();
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
void PETScSolver<Matrix_type>::factorizeMatrix()
{

  timer jactime; // start timing Jacobian evaluation

  SANS_ASSERT( continuousmap_.size() == 1); //for now...

  // update the matrix
  *A_ = 0;
  f_.jacobian(*A_);

  if (timing_ && continuousmap_[0].comm->rank() == 0) std::cout << "Jacobian time : " << jactime.elapsed() << " second(s)" << std::endl;

  // flatten the matrix for PETSc
  ScalarMatrix_CRS<PetscInt, PetscScalar> Ms(*A_);

#if 0
    std::string filename = "tmp/jac2.mtx";
    std::cout << "Writing Jacobian matrix to file: " << filename << "..." << std::endl;
    Ms.WriteMatrixMarketFile( filename );
#endif

  PetscScalar *a = Ms.Rx();

#ifdef SANS_MPI
  PetscInt m = Ms.m();
  PetscInt *row_ptr = Ms.Rp();
  PetscInt *col_ind = Ms.Ri();

  static const int M = StateVectorSize<SparseVectorView_type>::M;

  // convert the local column indexing to a global continuous index
  int nnz = Ms.nnz();
  for (int i = 0; i < nnz; i++)
  {
    if ( col_ind[i] < m )
      col_ind[i] += M*continuousmap_[0].nDOF_rank_offset;
    else
      col_ind[i] = col_ind[i] % M + M*continuousmap_[0].remoteGhostIndex[col_ind[i]/M - continuousmap_[0].nDOFpossessed];
  }

  PetscInt rstart, rend;
  PETSc_STATUS( MatGetOwnershipRange(A_petsc_, &rstart, &rend) );
  SANS_ASSERT( rstart == 0 + M*continuousmap_[0].nDOF_rank_offset );
  SANS_ASSERT( rend   == m + M*continuousmap_[0].nDOF_rank_offset );

  for (PetscInt row = 0; row < m; row++)
  {
    PetscInt rowGlobal = row + M*continuousmap_[0].nDOF_rank_offset;
    for (int k = row_ptr[row]; k < row_ptr[row+1]; k++ )
      PETSc_STATUS( MatSetValue(A_petsc_, rowGlobal, col_ind[k], a[k], INSERT_VALUES) );
  }

  // TODO: Not sure why this does not work...
  //for (PetscInt row = 0; row < m; row++)
  //  PETSc_STATUS( MatSetValuesRow(A_petsc_, row + M*continuousmap_[0].nDOF_rank_offset, a + row_ptr[row]) );

#else

  PetscScalar *a_petsc = pMs_->Rx();

  // PETSc uses the provided memory in pMs_
  int nnz = pMs_->nnz();
  for (int i = 0; i < nnz; i++)
    a_petsc[i] = a[i];

#endif

  factorize_petsc();
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
PetscErrorCode
PETScSolver<Matrix_type>::MatGetOrdering_MDF(Mat mat, MatOrderingType type, IS *irow, IS *icol)
{
  PetscFunctionBegin;

  SANS_ASSERT( strcmp(type, "mdf") == 0 );

  //Check that mat is a sequential matrix (since global MPIAIJ matrices forward the calls to each local matrix)
  MatType mattype;
  PETSc_STATUS( MatGetType(mat, &mattype) );
  SANS_ASSERT( strcmp(mattype, MATSEQAIJ) == 0 );

  PetscInt  nrow_compressed;
  PetscBool done;
  PETSc_STATUS( MatGetRowIJ(mat,0,PETSC_FALSE,PETSC_TRUE,&nrow_compressed,NULL,NULL,&done) );
  PETSc_STATUS( MatRestoreRowIJ(mat,0,PETSC_FALSE,PETSC_TRUE,NULL,NULL,NULL,&done) );

  PetscInt nrow = -1, ncol = -1;
  PETSc_STATUS( MatGetLocalSize(mat, &nrow, &ncol) );
  SANS_ASSERT(nrow > 0);
  SANS_ASSERT(nrow == ncol); //MPIAIJ matrices only forward the on-diagonal portion of local matrices, so should be square

  if (done && nrow != nrow_compressed)
    SANS_DEVELOPER_EXCEPTION("MDF ordering is not implemented for a symbolic factorization \"compressed\" due to i-nodes or block storage");

  PetscInt nCols;
  const PetscInt *cols;
  const PetscScalar *vals;

//  MPI_Comm comm;
//  PetscInt rank = 0;
//  PETSc_STATUS( PetscObjectGetComm((PetscObject)mat,&comm) );
//  MPI_Comm_rank(comm, &rank);

  //-----------------------------------------------
  //Create nonzero pattern of scalar matrix C

  SparseNonZeroPattern<Real> nz(nrow, ncol);

  for (int i = 0; i < nrow; i++)
  {
    PETSc_STATUS( MatGetRow(mat, i, &nCols, &cols, NULL) );

    for (int j = 0; j < nCols; j++)
      nz.add(i, cols[j]);

    PETSc_STATUS( MatRestoreRow(mat, i, &nCols, &cols, NULL) );
  }

  //-----------------------------------------------
  //Create C matrix : C(i,j) = | A(i,i)^-1 A(i,j) |

  SparseMatrix_CRS<Real> C(nz);
  C = 0;

  for (int i = 0; i < nrow; i++)
  {
    PETSc_STATUS( MatGetRow(mat, i, &nCols, &cols, &vals) );
    PetscScalar diag = 0;

    //Save off diagonal values
    for (int j = 0; j < nCols; j++)
    {
      if (i == cols[j])
      {
        diag = vals[j];
//        std::cout << "diag" << i << ": " << diag[i] << std::endl;
        break;
      }
    }

    //Compute Cij for off-diagonal entries
    for (int j = 0; j < nCols; j++)
    {
//      std::cout << "rank" << rank << ": (" << i << ", " << cols[j] << ") : " << vals[j] << std::endl;

      if (i != cols[j])
      {
        if (diag == 0)
        {
          C.sparseRow(i,j) = std::numeric_limits<Real>::max();
        }
        else
        {
          C.sparseRow(i,j) = fabs( vals[j]/diag ); //C(i,j) = | A(i,i)^-1 A(i,j) |
        }
//        std::cout << C.sparseRow(i,j) << std::endl;
      }
    }

    PETSc_STATUS( MatRestoreRow(mat, i, &nCols, &cols, &vals) );
  }

  std::vector<PetscInt> ordering = computeOrdering_MDF(C);

#if 0
  {
    SparseNonZeroPattern<Real> pivot_nz(nrow, ncol);

    for (int i = 0; i < nrow; i++)
      pivot_nz.add(i, ordering[i]);

    std::string filename = "tmp/pivot.mtx";
    std::cout << "Writing Pivot matrix to file: " << filename << "..." << std::endl;
    WriteMatrixMarketFile( pivot_nz, filename );
  }
#endif

  std::set<PetscInt> indexset;
  indexset.insert(ordering.begin(), ordering.end());

  SANS_ASSERT_MSG((PetscInt) indexset.size() == nrow,
                  "Index mapping is not unique! Indexset.size = %d, nrow = %d", indexset.size(), nrow );

  PETSc_STATUS( ISCreateGeneral(PETSC_COMM_SELF, nrow, ordering.data(), PETSC_COPY_VALUES, irow) );
  PETSc_STATUS( ISCreateGeneral(PETSC_COMM_SELF, nrow, ordering.data(), PETSC_COPY_VALUES, icol) );
  PetscFunctionReturn(0);

  // Taken from MatGetOrdering_Natural in src/mat/order/sorder.c
#if 0
  PetscErrorCode ierr;
  PetscInt       n,i,*ii;
  PetscBool      done;
  MPI_Comm       comm;

  PetscFunctionBegin;
  ierr = PetscObjectGetComm((PetscObject)mat,&comm);CHKERRQ(ierr);
  ierr = MatGetRowIJ(mat,0,PETSC_FALSE,PETSC_TRUE,&n,NULL,NULL,&done);CHKERRQ(ierr);
  ierr = MatRestoreRowIJ(mat,0,PETSC_FALSE,PETSC_TRUE,NULL,NULL,NULL,&done);CHKERRQ(ierr);
  if (done)
  { /* matrix may be "compressed" in symbolic factorization, due to i-nodes or block storage */
    /*
      We actually create general index sets because this avoids mallocs to
      to obtain the indices in the MatSolve() routines.
      ierr = ISCreateStride(PETSC_COMM_SELF,n,0,1,irow);CHKERRQ(ierr);
      ierr = ISCreateStride(PETSC_COMM_SELF,n,0,1,icol);CHKERRQ(ierr);
    */
    ierr = PetscMalloc1(n,&ii);CHKERRQ(ierr);
    for (i=0; i<n; i++) ii[i] = i;
    ierr = ISCreateGeneral(PETSC_COMM_SELF,n,ii,PETSC_COPY_VALUES,irow);CHKERRQ(ierr);
    ierr = ISCreateGeneral(PETSC_COMM_SELF,n,ii,PETSC_OWN_POINTER,icol);CHKERRQ(ierr);
  }
  else
  {
    PetscInt start,end;

    ierr = MatGetOwnershipRange(mat,&start,&end);CHKERRQ(ierr);
    ierr = ISCreateStride(comm,end-start,start,1,irow);CHKERRQ(ierr);
    ierr = ISCreateStride(comm,end-start,start,1,icol);CHKERRQ(ierr);
  }
  ierr = ISSetIdentity(*irow);CHKERRQ(ierr);
  ierr = ISSetIdentity(*icol);CHKERRQ(ierr);

  ISView(*irow, PETSC_VIEWER_STDOUT_SELF);
  ISView(*icol, PETSC_VIEWER_STDOUT_SELF);

  PetscFunctionReturn(0);
#endif

}

} //namespace SLA
} //namespace numpack 
