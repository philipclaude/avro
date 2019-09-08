// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "tools/timer.h"

#include "numpack/SparseLinAlg/ScalarVector.h"

#define PETSc_INSTANTIATE
#include "PETScSolver_impl.h"
#include "PETScSolver_factorize.h"

#define MPI_COMMUNICATOR_IN_CPP
#include "MPI/communicator.h"

#include <petscmat.h>

namespace numpack 
{
namespace SLA
{
#if 0
//-----------------------------------------------------------------------------
template< class Matrix_type >
void PETScSolver<Matrix_type>::syncGhosts( SparseVectorView_type& x ) const
{
#ifdef SANS_MPI
  Vec x_petsc;

  SANS_ASSERT( A_->n() == x.m() );

  PetscInt npossessed = continuousmap_.nDOFpossessed;
  PetscInt nghost = continuousmap_.remoteGhostIndex.size();
  const PetscInt *ghosts = continuousmap_.remoteGhostIndex.data();

  SANS_ASSERT( npossessed + nghost == x.m() );

  // Create the ghosted solution vector using the memory of x
  PETSc_STATUS( VecCreateGhostWithArray(*continuousmap_.comm, npossessed, PETSC_DECIDE, nghost, ghosts, &x[0], &x_petsc) );

  // update ghost information
  PETSc_STATUS( VecGhostUpdateBegin(x_petsc, INSERT_VALUES, SCATTER_FORWARD) );
  PETSc_STATUS( VecGhostUpdateEnd(x_petsc, INSERT_VALUES, SCATTER_FORWARD) );

  // release the memory from the PETSc vectors
  PETSc_STATUS( VecDestroy(&x_petsc) );
#endif
}
#endif

//-----------------------------------------------------------------------------
template< class Matrix_type >
LinearSolveStatus PETScSolver<Matrix_type>::
backsolve( const bool transpose, const SparseVectorView_type& b, SparseVectorView_type& x ) const
{
  Vec x_petsc, b_petsc;

  SANS_ASSERT( A_->m() == b.m() );
  SANS_ASSERT( A_->n() == x.m() );
  SANS_ASSERT( continuousmap_.size() == 1 );

#ifdef SANS_MPI

  PetscInt npossessed = continuousmap_[0].nDOFpossessed;
  PetscInt nghost = continuousmap_[0].remoteGhostIndex.size();
  const PetscInt *ghosts = continuousmap_[0].remoteGhostIndex.data();

  SANS_ASSERT_MSG( npossessed + nghost == x.m(), "npossessed=%d, nghost=%d, x.m()=%d", npossessed, nghost, x.m() );
  SANS_ASSERT( npossessed == b.m() );

  // Create the ghosted solution vector using the memory of x
  PETSc_STATUS( VecCreateGhostWithArray(*continuousmap_[0].comm, npossessed, PETSC_DECIDE, nghost, ghosts, &x[0], &x_petsc) );

  // Create the petsc residual vector using the memory of b
  PETSc_STATUS( VecCreateMPIWithArray(*continuousmap_[0].comm, 1, b.m(), PETSC_DECIDE, &b[0], &b_petsc) );

#else

  // Create a petsc solution vector using the memory of x
  PETSc_STATUS( VecCreateSeqWithArray(PETSC_COMM_SELF, 1, x.m(), &x[0], &x_petsc) );

  // Create the petsc residual vector using the memory of b
  PETSc_STATUS( VecCreateSeqWithArray(PETSC_COMM_SELF, 1, b.m(), &b[0], &b_petsc) );

#endif

  timer solvetime; // start timing solve

  // solve using petsc
  if (transpose)
  {
    PETSc_STATUS( KSPSolveTranspose(ksp_, b_petsc, x_petsc) );
  }
  else
  {
    PETSc_STATUS( KSPSolve(ksp_, b_petsc, x_petsc) );
  }

  if (timing_ && continuousmap_[0].comm->rank() == 0) std::cout << "PETSc solve time : " << solvetime.elapsed() << " second(s)" << std::endl;

  // update ghost information
  PETSc_STATUS( VecGhostUpdateBegin(x_petsc, INSERT_VALUES, SCATTER_FORWARD) );
  PETSc_STATUS( VecGhostUpdateEnd(x_petsc, INSERT_VALUES, SCATTER_FORWARD) );

  //Print debug info if needed
  getDebugInfo();

  // release the memory from the PETSc vectors
  PETSc_STATUS( VecDestroy(&b_petsc) );
  PETSc_STATUS( VecDestroy(&x_petsc) );

  KSPConvergedReason reason;
  PETSc_STATUS( KSPGetConvergedReason(ksp_, &reason) );

  return LinearSolveStatus(reason >= 0);
}

template class PETScSolver< SparseMatrix_CRS<double> >;

} //namespace SLA
} //namespace numpack 
