// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include <vector>
#include <algorithm> // std::sort
#include <set>

#include "tools/timer.h"
#include "tools/output_std_vector.h"

#define PETSc_INSTANTIATE
#include "PETScSolver_impl.h"
//#include "PETScSolver_factorize.h"
//#include "PETScSolver_backsolve.h"
#include "PETSc_VectorSize_impl.h"
#include "MDFOrdering.h"

#include "LinearAlgebra/DenseLinAlg/InverseLUP.h"
#include "LinearAlgebra/DenseLinAlg/DynamicSize/MatrixD_Norm.h"

#define SCALARMATRIX_CRS_INSTANTIATE
#include "LinearAlgebra/SparseLinAlg/ScalarMatrix_CRS_impl.h"

#define SCALARVECTOR_INSTANTIATE
#include "LinearAlgebra/SparseLinAlg/ScalarVector_impl.h"
#include "LinearAlgebra/SparseLinAlg/WritePlainVector.h"

#include "LinearAlgebra/SparseLinAlg/WriteMatrixMarketFile.h"
namespace SANS
{
namespace SLA
{


//-----------------------------------------------------------------------------
template< class Matrix_type >
LinearSolveStatus PETScSolver<Matrix_type>::solve(SparseVectorView_type& b, SparseVectorView_type& x)
{

  if ( staticCondensed_ )
  {
    factorize(b, transpose_);

    SANS_ASSERT( x.m() == 3 && b.m() == 3 ); //hardcode to static condensation from 3x3 system
    SparseVectorView_type bcondensed = b.sub(1,2);
    SparseVectorView_type xcondensed = x.sub(1,2);

    LinearSolveStatus status = backsolve(bcondensed, xcondensed);

    f_.completeUpdate(b, xcondensed, x);

    return status;
  }
  else
  {
    factorize();
    return backsolve( b, x );
  }

}

//-----------------------------------------------------------------------------
template< class Matrix_type >
void PETScSolver<Matrix_type>::init()
{
  if (memory_) { PETSc_STATUS( PetscMemorySetGetMaximumUsage() ); } //log memory consumption

  SystemNonZeroPattern nz(f_.matrixSize());

  if (staticCondensed_)
  {
    //dummy residual vector
    SparseVector_type dummyb(f_.vectorEqSize());
    dummyb = 0;

    f_.jacobian(dummyb, nz, transpose_);
  }
  else
  {
    f_.jacobian(nz);
  }

  // construct the matrix
  this->A_ = new Matrix_type(nz);

  SANS_ASSERT_MSG( (int)continuousmap_.size() == this->A_->n(),
                   "continuousmap_.size() = %d, this->A_->n() = %d",
                   continuousmap_.size(), this->A_->n() );
  SANS_ASSERT( this->A_->m() == this->A_->n() );

  DLA::VectorD< int > row_m(this->A_->m());
  DLA::VectorD< int > col_n(this->A_->n());

  row_m = 0;
  col_n = 0;

  PetscInt m = 0, n = 0, nnz = 0;

  for (int i = 0; i < this->A_->m(); i++)
  {
    for (int j = 0; j < this->A_->n(); j++)
    {
      row_m[i] = (*this->A_)(i,j).m();
      col_n[j] = (*this->A_)(i,j).n();

      nnz += (*this->A_)(i,j).getNumNonZero();
    }

    //Number of rows in scalar matrix
    m += row_m[i];
    n += col_n[i];
  }

  SANS_ASSERT( m <= n );  // Can only solve square matrices (but the matrix may contain off processor ghost columns)

  static const int M = StateVectorSize<SparseVectorView_type>::M;

  // Create the PETSc matrix
  std::vector<PetscInt> row_ind_petsc( m+1 );
  std::vector<PetscInt> col_ind_petsc( nnz );

  row_ind_petsc[0] = 0;
  int kk = 0;
  int row = 0;
  for (int i = 0; i < (*this->A_).m(); i++)
  {
    for ( int si = 0; si < row_m[i]; si++ )
    {
      row++;
      row_ind_petsc[row] = row_ind_petsc[row-1];

      int col = 0;
      for (int j = 0; j < (*this->A_).n(); j++)
      {
        if ( j > 0 ) col += continuousmap_[j-1].nDOFpossessed;

        typename Matrix_type::node_type& A = (*this->A_)(i,j);

        if (A.getNumNonZero() == 0) continue;

        PetscInt *col_ind = A.get_col_ind();
        PetscInt *row_ptr = A.get_row_ptr();

        row_ind_petsc[row] += row_ptr[si+1] - row_ptr[si];

        for ( int sk = row_ptr[si]; sk < row_ptr[si+1]; ++sk )
        {
          if ( col_ind[sk] < continuousmap_[j].nDOFpossessed )
            col_ind_petsc[kk] = col + col_ind[sk] + continuousmap_[j].nDOF_rank_offset;
          else
            col_ind_petsc[kk] = continuousmap_[j].remoteGhostIndex[col_ind[sk] - continuousmap_[j].nDOFpossessed];

          kk++;
        }
      }
    }
  }
  SANS_ASSERT( nnz == kk );

#if 0
  std::cout << "global -----" << std::endl;

  for (PetscInt row = 0; row < m; row++)
  {
    for (PetscInt sk = row_ind_petsc[row]; sk < row_ind_petsc[row+1]; sk++)
     std::cout << col_ind_petsc[sk] << ", ";
    std::cout << std::endl;
  }

  std::cout << "----------" << std::endl;
#endif

#ifdef SANS_MPI

  // PETSc requires the columns to be sorted for block matrices (it's a bug...)
  if (nnz > 0) // suppress analyzer warnings
    for (PetscInt row = 0; row < m; row++)
      std::sort(col_ind_petsc.begin() + row_ind_petsc[row], col_ind_petsc.begin() + row_ind_petsc[row+1]);

  // So there seems to be a bug in this function, below is the corrected guts
  //PETSc_STATUS( MatCreateMPIBAIJWithArrays(*continuousmap_.comm, M, M*m, M*m, PETSC_DETERMINE, PETSC_DETERMINE,
  //                                         row_ptr, col_ind_petsc.data(), nullptr, &A_petsc_) );

  // Correct guts for MatCreateMPIBAIJWithArrays
  PETSc_STATUS( MatCreate(*continuousmap_[0].comm, &A_petsc_) );
  PETSc_STATUS( MatSetSizes(A_petsc_, M*m, M*m, PETSC_DETERMINE, PETSC_DETERMINE) );
  PETSc_STATUS( MatSetType(A_petsc_, MATMPIBAIJ) );
  PETSc_STATUS( MatSetBlockSize(A_petsc_, M) );
  PETSc_STATUS( MatSetUp(A_petsc_) );
  PETSc_STATUS( MatMPIBAIJSetPreallocationCSR(A_petsc_, M, row_ind_petsc.data(), col_ind_petsc.data(), nullptr) );

#else

  PETSc_STATUS( MatCreate(*continuousmap_[0].comm, &A_petsc_) );
  PETSc_STATUS( MatSetSizes(A_petsc_, M*m, M*n, PETSC_DETERMINE, PETSC_DETERMINE) );
  PETSc_STATUS( MatSetType(A_petsc_,MATSEQBAIJ) );
  PETSc_STATUS( MatSeqBAIJSetPreallocationCSR(A_petsc_, M, row_ind_petsc.data(), col_ind_petsc.data(), nullptr) );
#endif

  init_petsc();
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
void PETScSolver<Matrix_type>::
factorize(SparseVectorView_type& b, bool transpose)
{
  SANS_ASSERT(staticCondensed_ == true);
  timer jactime; // start timing Jacobian evaluation

  // update the matrix
  *A_ = 0;

//  const std::vector<std::vector<Real>> nrmRsd1 = f_.residualNorm(b);
//
//  std::string file1 = "tmp/rsd1_rank";
//  file1 += std::to_string(continuousmap_[0].comm->rank());
//  file1 += ".dat";
//  WritePlainVector( b, file1 );

  f_.jacobian(b, *A_, transpose);

//  const std::vector<std::vector<Real>> nrmRsd2 = f_.residualNorm(b);
//
//  std::string file2 = "tmp/rsd2_rank";
//  file2 += std::to_string(continuousmap_[0].comm->rank());
//  file2 += ".dat";
//  WritePlainVector( b, file2 );


  if (timing_ && continuousmap_[0].comm->rank() == 0) std::cout << "Jacobian time : " << jactime.elapsed() << " second(s)" << std::endl;

  factorizeMatrix();
}


//-----------------------------------------------------------------------------
template< class Matrix_type >
void PETScSolver<Matrix_type>::factorize()
{

  SANS_ASSERT(staticCondensed_ == false);
  timer jactime; // start timing Jacobian evaluation

  // update the matrix
  *A_ = 0;

  f_.jacobian(*A_);
  if (timing_ && continuousmap_[0].comm->rank() == 0) std::cout << "Jacobian time : " << jactime.elapsed() << " second(s)" << std::endl;

  factorizeMatrix();
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
void PETScSolver<Matrix_type>::factorizeMatrix()
{

  DLA::VectorD< int > row_m(this->A_->m());
  DLA::VectorD< int > col_n(this->A_->n());

  row_m = 0;
  col_n = 0;

  PetscInt m = 0, n = 0, nnz = 0;

  for (int i = 0; i < this->A_->m(); i++)
  {
    for (int j = 0; j < this->A_->n(); j++)
    {
      row_m[i] = (*this->A_)(i,j).m();
      col_n[j] = (*this->A_)(i,j).n();

      nnz += (*this->A_)(i,j).getNumNonZero();
    }

    //Number of rows in scalar matrix
    m += row_m[i];
    n += col_n[i];
  }

  SANS_ASSERT( m <= n );  // Can only solve square matrices (but the matrix may contain off processor ghost columns)

  // Create the PETSc matrix
  std::vector<PetscInt> row_ind_petsc( m+1 );
  std::vector<PetscInt> col_ind_petsc( nnz );
  std::vector<typename Matrix_type::node_type::Ttype> a( nnz );

  row_ind_petsc[0] = 0;
  int kk = 0;
  int row = 0;
  for (int i = 0; i < (*this->A_).m(); i++)
  {
    for ( int si = 0; si < row_m[i]; si++ )
    {
      row++;
      row_ind_petsc[row] = row_ind_petsc[row-1];

      int col = 0;
      for (int j = 0; j < (*this->A_).n(); j++)
      {
        if ( j > 0 ) col += continuousmap_[j-1].nDOFpossessed;

        typename Matrix_type::node_type& A = (*this->A_)(i,j);

        if (A.getNumNonZero() == 0) continue;

        PetscInt *col_ind = A.get_col_ind();
        PetscInt *row_ptr = A.get_row_ptr();

        row_ind_petsc[row] += row_ptr[si+1] - row_ptr[si];

        for ( int sk = row_ptr[si]; sk < row_ptr[si+1]; ++sk )
        {
          if ( col_ind[sk] < continuousmap_[j].nDOFpossessed )
            col_ind_petsc[kk] = col + col_ind[sk] + continuousmap_[j].nDOF_rank_offset;
          else
            col_ind_petsc[kk] = continuousmap_[j].remoteGhostIndex[col_ind[sk] - continuousmap_[j].nDOFpossessed];

          a[kk] = A[sk];
          kk++;
        }
      }
    }
  }

#if 0
  std::cout << "global -----" << std::endl;

  for (PetscInt row = 0; row < m; row++)
  {
    for (PetscInt sk = row_ind_petsc[row]; sk < row_ind_petsc[row+1]; sk++)
     std::cout << col_ind_petsc[sk] << ", ";
    std::cout << std::endl;
  }

  std::cout << "----------" << std::endl;
#endif

#ifdef SANS_MPI

  static const int M = StateVectorSize<SparseVectorView_type>::M;

  PetscInt rstart, rend;
  PETSc_STATUS( MatGetOwnershipRange(A_petsc_, &rstart, &rend) );
  SANS_ASSERT( rstart == M*(0 + continuousmap_[0].nDOF_rank_offset) );
  SANS_ASSERT( rend   == M*(m + continuousmap_[0].nDOF_rank_offset) );

  // copy over the memory
  if (nnz > 0) // suppress analyzer warning
    for (PetscInt row = 0; row < m; row++)
    {
      PetscInt rowGlobal = row + continuousmap_[0].nDOF_rank_offset;

      for (int k = row_ind_petsc[row]; k < row_ind_petsc[row+1]; k++ )
        PETSc_STATUS( MatSetValuesBlocked(A_petsc_, 1, &rowGlobal, 1, &col_ind_petsc[k], &DLA::index(a[k],0,0), INSERT_VALUES) );
    }

#else

  // copy over the memory
  if (nnz > 0) // suppress analyzer warning
    for (PetscInt row = 0; row < m; row++)
      for (int k = row_ind_petsc[row]; k < row_ind_petsc[row+1]; k++ )
        PETSc_STATUS( MatSetValuesBlocked(A_petsc_, 1, &row, 1, &col_ind_petsc[k], &DLA::index(a[k],0,0), INSERT_VALUES) );

#endif

  factorize_petsc();
}

#if 0
//-----------------------------------------------------------------------------
template< class Matrix_type >
void PETScSolver<Matrix_type>::syncGhosts( SparseVectorView_type& x ) const
{
#ifdef SANS_MPI
  Vec x_petsc;

  SANS_ASSERT( A_->n() == x.m() );

  // TODO: Still need to implement for MatrixD larger than 1x1
  for (int i = 1; i < x.m(); i++)
    SANS_ASSERT(x[i].m() == 0);

  static const int M = StateVectorSize<SparseVectorView_type>::M;

  PetscInt n = x[0].m()*M;

  std::vector<PetscScalar> xx( n );

  for ( int i = 0; i < x[0].m(); i++ )
    for ( int vi = 0; vi < M; vi++ )
      xx[M*i+vi] = x[0][i][vi];

  PetscInt npossessed = continuousmap_.nDOFpossessed;
  PetscInt nghost = continuousmap_.remoteGhostIndex.size();
  std::vector<PetscInt> ghosts(nghost);
  for (std::size_t i = 0; i < continuousmap_.remoteGhostIndex.size(); i++)
    ghosts[i] = continuousmap_.remoteGhostIndex[i];

  SANS_ASSERT( npossessed + nghost == x[0].m() );

  // Create the ghosted solution vector using the memory of x
  PETSc_STATUS( VecCreateGhostBlockWithArray(*continuousmap_.comm, M, M*npossessed, PETSC_DECIDE, nghost, ghosts.data(), xx.data(), &x_petsc) );

  // update ghost information
  PETSc_STATUS( VecGhostUpdateBegin(x_petsc, INSERT_VALUES, SCATTER_FORWARD) );
  PETSc_STATUS( VecGhostUpdateEnd(x_petsc, INSERT_VALUES, SCATTER_FORWARD) );

  // copy the information back from the scalar vector
  for ( int i = 0; i < x[0].m(); i++ )
    for ( int vi = 0; vi < M; vi++ )
      x[0][i][vi] = xx[M*i+vi];

  // release the memory from the PETSc vectors
  PETSc_STATUS( VecDestroy(&x_petsc) );
#endif
}
#endif

//-----------------------------------------------------------------------------
template< class Matrix_type >
LinearSolveStatus PETScSolver<Matrix_type>::
backsolve(const bool transpose, const SparseVectorView_type& b, SparseVectorView_type& x ) const
{
  Vec x_petsc, b_petsc;

  SANS_ASSERT_MSG( A_->m() == b.m(), "A_->m() = %d, b.m() = %d", A_->m(), b.m() );
  SANS_ASSERT_MSG( A_->n() == x.m(), "A_->m() = %d, x.m() = %d", A_->m(), x.m() );

  static const int M = StateVectorSize<SparseVectorView_type>::M;

  PetscInt m = 0;
  PetscInt n = 0;

  for (int i = 0; i < b.m(); i++)
    m += b[i].m();

  for (int j = 0; j < x.m(); j++)
    n += x[j].m();

  std::vector<PetscScalar> bb( M*m );
  std::vector<PetscScalar> xx( M*n );

  int ii = 0;
  for ( int i = 0; i < x.m(); i++ )
    for ( int si = 0; si < b[i].m(); si++ ) // must use b.m here go get possessed DOFs
      for ( int bi = 0; bi < M; bi++ )
        xx[ii++] = DLA::index(x[i][si],bi);

  // set ghost DOFs
  for ( int i = 0; i < x.m(); i++ )
    for ( int si = b[i].m(); si < x[i].m(); si++ )
      for ( int bi = 0; bi < M; bi++ )
        xx[ii++] = DLA::index(x[i][si],bi);

  ii = 0;
  for ( int i = 0; i < b.m(); i++ )
    for ( int si = 0; si < b[i].m(); si++ )
      for ( int bi = 0; bi < M; bi++ )
        bb[ii++] = DLA::index(b[i][si],bi);

#ifdef SANS_MPI


  PetscInt npossessed = 0;
  PetscInt nghost = 0;
  for (std::size_t iv = 0; iv < continuousmap_.size(); iv++)
  {
    npossessed += continuousmap_[iv].nDOFpossessed;
    nghost     += continuousmap_[iv].remoteGhostIndex.size();
  }

  ii = 0;
  std::vector<PetscInt> ghosts(nghost);
  for (std::size_t iv = 0; iv < continuousmap_.size(); iv++)
    for (std::size_t i = 0; i < continuousmap_[iv].remoteGhostIndex.size(); i++)
    {
      ghosts[ii++] = continuousmap_[iv].remoteGhostIndex[i];
    }

  SANS_ASSERT_MSG( npossessed + nghost == n,
                   "npossessed=%d, nghost=%d, n=%d, rank=%d",
                   npossessed, nghost, n, continuousmap_[0].comm->rank() );
  SANS_ASSERT( npossessed == m );

  // Create the ghosted solution vector using the memory of x
  PETSc_STATUS( VecCreateGhostBlockWithArray(*continuousmap_[0].comm, M, M*npossessed, PETSC_DECIDE, nghost, ghosts.data(), xx.data(), &x_petsc) );

  // Create the petsc residual vector using the memory of b
  PETSc_STATUS( VecCreateMPIWithArray(*continuousmap_[0].comm, M, M*m, PETSC_DECIDE, bb.data(), &b_petsc) );

#else

  // Create a petsc solution vector using the memory of x
  PETSc_STATUS( VecCreateSeqWithArray(PETSC_COMM_SELF, M, M*n, xx.data(), &x_petsc) );

  // Create the petsc residual vector using the memory of b
  PETSc_STATUS( VecCreateSeqWithArray(PETSC_COMM_SELF, M, M*m, bb.data(), &b_petsc) );

#endif

  timer solvetime; // start timing solve

#if 0 // Show properties about how the KSP is set
  PetscViewer viewer;
  PETSc_STATUS(PetscViewerCreate(*continuousmap_.comm, &viewer));
  PETSc_STATUS(PetscViewerSetType(viewer,PETSCVIEWERASCII));
  PETSc_STATUS(KSPView(ksp_, viewer));
  PETSc_STATUS(PetscViewerDestroy(&viewer));

#endif

  // solve using petsc
  if (transpose)
  {
#if 0 // Dump matrix and rhs for debugging
    PetscViewer matrixfile, rhsfile;
    PETSc_STATUS(PetscViewerBinaryOpen(*continuousmap_[0].comm, "tmp/matrix.dat", FILE_MODE_WRITE, &matrixfile ));

    PETSc_STATUS(MatView(A_petsc_,matrixfile));
    PETSc_STATUS(PetscViewerDestroy(&matrixfile));

    PETSc_STATUS(PetscViewerBinaryOpen(*continuousmap_[0].comm, "tmp/rhs.dat", FILE_MODE_WRITE, &rhsfile ));
    PETSc_STATUS(VecView(b_petsc,rhsfile));
    PETSc_STATUS(PetscViewerDestroy(&rhsfile));
#endif
    PETSc_STATUS( KSPSolveTranspose(ksp_, b_petsc, x_petsc) );
  }
  else
  {
    PETSc_STATUS( KSPSolve(ksp_, b_petsc, x_petsc) );
  }

  if (timing_ && continuousmap_[0].comm->rank() == 0) std::cout << "PETSc solve time : " << solvetime.elapsed() << " second(s)" << std::endl;

  if (memory_)
  {
    PetscLogDouble mem = 0;
    PETSc_STATUS( PetscMemoryGetCurrentUsage(&mem) );

#ifdef SANS_MPI
    mem = boost::mpi::all_reduce(*continuousmap_[0].comm, mem, std::plus<PetscLogDouble>());
#endif

    if (continuousmap_[0].comm->rank() == 0) std::cout << "PETSc Memory Usage : " << mem/(1024*1024) << " MB" << std::endl;

  }

  // update ghost information
  PETSc_STATUS( VecGhostUpdateBegin(x_petsc, INSERT_VALUES, SCATTER_FORWARD) );
  PETSc_STATUS( VecGhostUpdateEnd(x_petsc, INSERT_VALUES, SCATTER_FORWARD) );

  //Print debug info if needed
  getDebugInfo();

  // copy the information back from the scalar vector
  ii = 0;
  for ( int i = 0; i < x.m(); i++ )
    for ( int si = 0; si < b[i].m(); si++ ) // must be b.m to get possessed DOFs
      for ( int bi = 0; bi < M; bi++ )
        DLA::index(x[i][si],bi) = xx[ii++];

  for ( int i = 0; i < x.m(); i++ )
    for ( int si = b[i].m(); si < x[i].m(); si++ )
      for ( int bi = 0; bi < M; bi++ )
        DLA::index(x[i][si],bi) = xx[ii++];

  // release the memory from the PETSc vectors
  PETSc_STATUS( VecDestroy(&b_petsc) );
  PETSc_STATUS( VecDestroy(&x_petsc) );

  KSPConvergedReason reason;
  PETSc_STATUS( KSPGetConvergedReason(ksp_, &reason) );

  return LinearSolveStatus(reason >= 0);
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
PetscErrorCode
PETScSolver<Matrix_type>::MatGetOrdering_MDF(Mat mat, MatOrderingType type, IS *irow, IS *icol)
{
  PetscFunctionBegin;

  SANS_ASSERT( strcmp(type, "mdf") == 0 );

  //Check that mat is a sequential block matrix (since global MPIBAIJ matrices forward the calls to each local matrix)
  MatType mattype;
  PETSc_STATUS( MatGetType(mat, &mattype) );
  SANS_ASSERT( strcmp(mattype, MATSEQBAIJ) == 0 );

//  MPI_Comm comm;
//  PetscInt rank = 0;
//  PETSc_STATUS( PetscObjectGetComm((PetscObject)mat,&comm) );
//  MPI_Comm_rank(comm, &rank);

  PetscInt  nrow_compressed;
  PetscBool done;
  PETSc_STATUS( MatGetRowIJ(mat,0,PETSC_FALSE,PETSC_TRUE,&nrow_compressed,NULL,NULL,&done) );
  PETSc_STATUS( MatRestoreRowIJ(mat,0,PETSC_FALSE,PETSC_TRUE,NULL,NULL,NULL,&done) );

  PetscInt nrow = -1, ncol = -1;
  PETSc_STATUS( MatGetLocalSize(mat, &nrow, &ncol) );
  SANS_ASSERT(nrow > 0);
  SANS_ASSERT(nrow == ncol); //only called with on-diagonal part of local matrices, so should be square

  PetscContainer container;
  PETSc_STATUS( PetscObjectQuery((PetscObject) mat, "MDFouterblocksize", (PetscObject*) &container) );
  void * vptr;
  PETSc_STATUS( PetscContainerGetPointer(container, &vptr) );
  const int outer_block_size = *((int*) vptr);
  SANS_ASSERT( outer_block_size > 0 );
//  std::cout << "MDFOuterBlockSize: " << outer_block_size << std::endl;

  static const int M = StateVectorSize<SparseVectorView_type>::M;
  const int block_size = M*outer_block_size;

  SANS_ASSERT(nrow % block_size == 0); //ensure that there is an integer number of blocks
  const int nBlock = nrow / block_size;

  if (done && nrow/M != nrow_compressed)
    SANS_DEVELOPER_EXCEPTION("MDF ordering is not implemented for a symbolic factorization \"compressed\" due to i-nodes or block storage");

  PetscInt nCols;
  const PetscInt *cols;
  const PetscScalar *vals;

  //-----------------------------------------------
  //Create nonzero pattern of scalar matrix C

  SparseNonZeroPattern<Real> nz(nBlock, nBlock);

  for (int ib = 0; ib < nBlock; ib++)
  {
    for (int m = 0; m < block_size; m++)
    {
      const int i = block_size*ib + m;
      PETSc_STATUS( MatGetRow(mat, i, &nCols, &cols, NULL) );

      for (int j = 0; j < nCols; j++)
      {
        const int jb = cols[j] / block_size; //get the block column index
//        std::cout << "i: " << i << ", j: " << cols[j] << ", jb: " << jb << std::endl;
        nz.add(ib, jb);
      }

      PETSc_STATUS( MatRestoreRow(mat, i, &nCols, &cols, NULL) );
    }
  }

  //-----------------------------------------------
  //Create C matrix : C(i,j) = | A(i,i)^-1 A(i,j) |

  SparseMatrix_CRS<Real> C(nz);
  SANS_ASSERT( C.getNumNonZero() > 0 );
  C = 0;

  DLA::MatrixD<Real> zero_block(block_size, block_size);
  zero_block = 0.0;

  for (int ib = 0; ib < nBlock; ib++)
  {
    //Create a vector of MatrixDs for each "potential" block in this row (some blocks may not be filled)
    std::vector<DLA::MatrixD<Real>> block_row(nBlock, zero_block);

    //Get data for each block in this row
    for (int m = 0; m < block_size; m++)
    {
      const int i = block_size*ib + m;
      PETSc_STATUS( MatGetRow(mat, i, &nCols, &cols, &vals) );

      for (int j = 0; j < nCols; j++)
      {
        const int n = cols[j] % block_size; //column offset inside block
        const int jb = cols[j] / block_size; //get the block column index
//        std::cout << "i: " << i << ", j: " << cols[j] << ", jb: " << jb << std::endl;
        block_row[jb](m,n) = vals[j];
      }

      PETSc_STATUS( MatRestoreRow(mat, i, &nCols, &cols, &vals) );
    }

    // factorize block_row[ib]
    auto block_row_ib_fac = DLA::InverseLUP::Factorize(block_row[ib]);

    for (int j = 0; j < C.rowNonZero(ib); j++)
    {
      const int jb = C.columnIndex(ib, j); //get block column index

      if (ib != jb)
      {
        block_row[jb] = block_row_ib_fac.backsolve(block_row[jb]); // A(i,i)^-1 A(i,j)
        C.sparseRow(ib,j) = normFrobenius(block_row[jb]);
//        std::cout << "ib: " << ib << ", jb: " << jb << ", C: " << C.sparseRow(ib,j) << std::endl;
//        std::cout << block_row[jb] << std::endl;
      }
    }
  }

  //Compute the ordering of the outer blocks
  std::vector<int> ordering_scalar = computeOrdering_MDF(C);

  const int nBlock_PETSc = nrow / M;
  std::vector<int> ordering(nBlock_PETSc); //length of ordering should be the number of PETSc blocks (not total rows)
  std::set<int> indexset;

//  std::cout << "Ordering:" << std::endl;
  for (int ib = 0; ib < nBlock; ib++)
    for (int m = 0; m < outer_block_size; m++)
    {
      int orig_block_index = outer_block_size*ib + m;
      int new_block_index = outer_block_size*ordering_scalar[ib] + m;
      SANS_ASSERT(new_block_index >= 0 && new_block_index < nBlock_PETSc);

      ordering[orig_block_index] = new_block_index;
      indexset.insert(ordering[orig_block_index]);
//      std::cout << ordering[orig_block_index] << std::endl;
    }

  SANS_ASSERT_MSG((int) indexset.size() == nBlock_PETSc,
                  "Index mapping is not unique! Indexset.size = %d, nBlock = %d", indexset.size(), nBlock_PETSc );

  PETSc_STATUS( ISCreateGeneral(PETSC_COMM_SELF, nBlock_PETSc, ordering.data(), PETSC_COPY_VALUES, irow) );
  PETSc_STATUS( ISCreateGeneral(PETSC_COMM_SELF, nBlock_PETSc, ordering.data(), PETSC_COPY_VALUES, icol) );
  PetscFunctionReturn(0);
}

//-----------------------------------------------------------------------------
template class PETScSolver< DLA::MatrixD< SparseMatrix_CRS< Real> > >;

template class PETScSolver< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixS<2,2,Real> > > >;
template class PETScSolver< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixS<3,3,Real> > > >;
template class PETScSolver< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixS<4,4,Real> > > >;
template class PETScSolver< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixS<5,5,Real> > > >;
template class PETScSolver< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixS<6,6,Real> > > >;
template class PETScSolver< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixS<7,7,Real> > > >;
template class PETScSolver< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixS<8,8,Real> > > >;

template ScalarMatrix_CRS<PetscInt>::ScalarMatrix_CRS( const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<2,2,Real> > >&);
template ScalarMatrix_CRS<PetscInt>::ScalarMatrix_CRS( const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<3,3,Real> > >&);
template ScalarMatrix_CRS<PetscInt>::ScalarMatrix_CRS( const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<4,4,Real> > >&);
template ScalarMatrix_CRS<PetscInt>::ScalarMatrix_CRS( const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<5,5,Real> > >&);
template ScalarMatrix_CRS<PetscInt>::ScalarMatrix_CRS( const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<6,6,Real> > >&);
template ScalarMatrix_CRS<PetscInt>::ScalarMatrix_CRS( const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<7,7,Real> > >&);
template ScalarMatrix_CRS<PetscInt>::ScalarMatrix_CRS( const DLA::MatrixDView< SparseMatrix_CRS< DLA::MatrixS<8,8,Real> > >&);

} //namespace SLA
} //namespace SANS
