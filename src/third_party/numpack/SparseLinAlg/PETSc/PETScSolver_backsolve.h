// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(PETSc_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "tools/timer.h"

#include "numpack/SparseLinAlg/ScalarVector.h"

#include "PETScSolver.h"
#include "PETSc_VectorSize_impl.h"

#define MPI_COMMUNICATOR_IN_CPP
#include "MPI/communicator.h"

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

  ScalarVector xx(x);

  static const int M = StateVectorSize<SparseVectorView_type>::M;

  PetscInt npossessed = M*continuousmap_.nDOFpossessed;
  PetscInt nghost = M*continuousmap_.remoteGhostIndex.size();
  std::vector<PetscInt> ghosts(nghost);
  for (std::size_t i = 0; i < continuousmap_.remoteGhostIndex.size(); i++)
  {
    int remoteGlobalIndex = continuousmap_.remoteGhostIndex[i];
    for (int k = 0; k < M; k++)
      ghosts[M*i + k] = M*remoteGlobalIndex + k;
  }

  SANS_ASSERT( npossessed + nghost == xx.m );

  // Create the ghosted solution vector using the memory of x
  PETSc_STATUS( VecCreateGhostWithArray(*continuousmap_.comm, npossessed, PETSC_DECIDE, nghost, ghosts.data(), xx.v, &x_petsc) );

  // update ghost information
  PETSc_STATUS( VecGhostUpdateBegin(x_petsc, INSERT_VALUES, SCATTER_FORWARD) );
  PETSc_STATUS( VecGhostUpdateEnd(x_petsc, INSERT_VALUES, SCATTER_FORWARD) );

  // copy the information back from the scalar vector
  xx.setTo(x);

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

  ScalarVector xx(x);
  ScalarVector bb(b);

#ifdef SANS_MPI

  static const int M = StateVectorSize<SparseVectorView_type>::M;

  PetscInt npossessed = 0;
  PetscInt nghost = 0;
  for (std::size_t iv = 0; iv < continuousmap_.size(); iv++)
  {
    npossessed += M*continuousmap_[iv].nDOFpossessed;
    nghost     += M*continuousmap_[iv].remoteGhostIndex.size();
  }

  int ii = 0;
  std::vector<PetscInt> ghosts(nghost);
  for (std::size_t iv = 0; iv < continuousmap_.size(); iv++)
    for (std::size_t i = 0; i < continuousmap_[iv].remoteGhostIndex.size(); i++)
    {
      int remoteGlobalIndex = continuousmap_[iv].remoteGhostIndex[i];
      for (int k = 0; k < M; k++)
        ghosts[ii++] = M*remoteGlobalIndex + k;
    }

  SANS_ASSERT( npossessed + nghost == xx.m );
  SANS_ASSERT( npossessed == bb.m );

  // Create the ghosted solution vector using the memory of x
  PETSc_STATUS( VecCreateGhostWithArray(*continuousmap_[0].comm, npossessed, PETSC_DECIDE, nghost, ghosts.data(), xx.v, &x_petsc) );

  // Create the petsc residual vector using the memory of b
  PETSc_STATUS( VecCreateMPIWithArray(*continuousmap_[0].comm, 1, bb.m, PETSC_DECIDE, bb.v, &b_petsc) );

#else

  // Create a petsc solution vector using the memory of x
  PETSc_STATUS( VecCreateSeqWithArray(PETSC_COMM_SELF, 1, xx.m, xx.v, &x_petsc) );

  // Create the petsc residual vector using the memory of b
  PETSc_STATUS( VecCreateSeqWithArray(PETSC_COMM_SELF, 1, bb.m, bb.v, &b_petsc) );

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

  // copy the information back from the scalar vector
  xx.setTo(x);

  // release the memory from the PETSc vectors
  PETSc_STATUS( VecDestroy(&b_petsc) );
  PETSc_STATUS( VecDestroy(&x_petsc) );

  KSPConvergedReason reason;
  PETSc_STATUS( KSPGetConvergedReason(ksp_, &reason) );

  return LinearSolveStatus(reason >= 0);
}

} // namespace SLA
} // namespace numpack 
