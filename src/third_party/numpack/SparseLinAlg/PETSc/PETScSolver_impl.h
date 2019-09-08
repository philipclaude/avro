// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(PETSc_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "PETScSolver.h"

#define MPI_COMMUNICATOR_IN_CPP
#include "MPI/communicator.h"

#ifdef SANS_MPI
#include <boost/mpi/collectives/all_reduce.hpp>
#endif

#include <fstream>
#include <limits>
#include <algorithm>

#include <petscoptions.h>
#include <boost/filesystem/operations.hpp>

namespace numpack 
{
namespace SLA
{

//-----------------------------------------------------------------------------
template< class Matrix_type >
PETScSolver<Matrix_type>::ResidualMonitor::
ResidualMonitor(std::string filename, bool verbose, PetscInt MaxIterations, int comm_rank_) :
  dumpResidualHistory(!filename.empty()), verbose(verbose), MaxIterations(MaxIterations), comm_rank(comm_rank_)
{
  if (dumpResidualHistory && comm_rank == 0)
  {
    // check the see if the file already exists before opening it
    bool writeHeader = !boost::filesystem::exists( filename );

    fhist.open(filename, std::ios_base::app);

    // first time opening the file, add the header
    if (writeHeader)
      fhist << "VARIABLES=Iteration, Residual" << std::endl;

    // add a new zone to distinguish between restarts
    fhist << "ZONE" << std::endl;
  }
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
PETScSolver<Matrix_type>::PETScSolver( PyDict d, AlgebraicEquationSetBase<Matrix_type>& f, LinearSystemSolve solve ) :
  Base_type(solve),
  params(PETScSolverParam::params),
  f_(f),
  continuousmap_(f.continuousGlobalMap()),
  pMs_(nullptr),
  A_petsc_(nullptr),
  ksp_(nullptr),
  RelativeTolerance_(d.get(params.RelativeTolerance)),
  AbsoluteTolerance_(d.get(params.AbsoluteTolerance)),
  DivergenceTolerance_(d.get(params.DivergenceTolerance)),
  MaxIterations_(d.get(params.MaxIterations)),
  GMRES_restart_(d.get(params.GMRES_Restart)),
  PreconditionerDict_(d.get(params.Preconditioner)),
  verbose_(d.get(params.Verbose)),
  timing_(d.get(params.Timing)),
  memory_(d.get(params.Memory)),
  computeSingularValues_(d.get(params.computeSingularValues)),
  printMatrixInfo_(d.get(params.printMatrixInfo)),
  matprefix_(d.get(params.PETScMatPrefix)),
  kspprefix_(d.get(params.PETScKSPPrefix)),
  rsdMon_(d.get(params.ResidualHistoryFile), verbose_, MaxIterations_,
          continuousmap_[0].comm != nullptr ? continuousmap_[0].comm->rank() : 0),
  solverType_(d.get(params.KSPSolver)),
  staticCondensed_(f.isStaticCondensed())
{
#ifdef SANS_MPI
  SANS_ASSERT( continuousmap_[0].comm != nullptr );
#endif

  // create a PETSc options database object
  std::string PETScOptions = d.get(params.PETScOptions);
  std::string PETScOptionsFile = d.get(params.PETScOptionsFile);

  PETSc_STATUS( PetscOptionsInsertString(NULL, PETScOptions.c_str()) );
  if (!PETScOptionsFile.empty())
  {
#ifdef SANS_MPI
    PETSc_STATUS( PetscOptionsInsertFile(*continuousmap_[0].comm, NULL, PETScOptionsFile.c_str(), PETSC_TRUE) );
#else
    PETSc_STATUS( PetscOptionsInsertFile(PETSC_COMM_SELF, NULL, PETScOptionsFile.c_str(), PETSC_TRUE) );
#endif
  }

  PETSc_STATUS( PetscContainerCreate(*continuousmap_[0].comm, const_cast<PetscContainer*>(&mdfcontainer_) ) );

  init();
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
void PETScSolver<Matrix_type>::init_petsc()
{
  PETSc_STATUS( PetscObjectSetOptionsPrefix((PetscObject)A_petsc_, matprefix_.c_str()) );
  PETSc_STATUS( MatSetFromOptions(A_petsc_) );

  // the matrix is done, and setup locally on each processor
  PETSc_STATUS( MatSetOption(A_petsc_, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE) );
  PETSc_STATUS( MatSetOption(A_petsc_, MAT_NO_OFF_PROC_ENTRIES  , PETSC_TRUE) );

  PETSc_STATUS( MatSetUp(A_petsc_) );

#if 0
  if (verbose_)
  {
    MatType type;
    PetscInt mout, nout;
    PETSc_STATUS( MatGetType(A_petsc_, &type) );
    PETSc_STATUS( MatGetSize(A_petsc_, &mout, &nout) );
    std::cout << "PETSc matrix type: " << type << ", size: " << mout << " x " << nout << std::endl;
  }
#endif

  // Create the linear solver
  PETSc_STATUS( KSPCreate(*continuousmap_[0].comm, &ksp_) );
  PETSc_STATUS( PetscObjectSetOptionsPrefix((PetscObject)ksp_, kspprefix_.c_str()) );

  // Set the matrix for the linear solver
  PETSc_STATUS( KSPSetOperators(ksp_, A_petsc_, A_petsc_) );

  // Set up linear solve parameters
  PETSc_STATUS( KSPSetTolerances(ksp_, RelativeTolerance_, AbsoluteTolerance_, DivergenceTolerance_, MaxIterations_) );

  if (rsdMon_.verbose || rsdMon_.dumpResidualHistory)
    PETSc_STATUS( KSPMonitorSet(ksp_, residualMonitor, &rsdMon_, nullptr) );

  if (computeSingularValues_)
    PETSc_STATUS( KSPSetComputeSingularValues(ksp_, PETSC_TRUE) );

  // Use current values in the x vector as initial guess (i.e. don't set x = 0 inside PETSc)
  PETSc_STATUS( KSPSetInitialGuessNonzero(ksp_,PETSC_TRUE) );

  //set any additional options
  PETSc_STATUS( KSPSetFromOptions(ksp_) );

  if (solverType_ == PETScSolverParam::params.KSPSolver.GMRES) {}
  else if (solverType_ == PETScSolverParam::params.KSPSolver.BICGStab)
  {
    PETSc_STATUS( KSPSetType(ksp_, KSPBCGSL) );
  }
  else if (solverType_ == PETScSolverParam::params.KSPSolver.DGMRES)
  {
    PETSc_STATUS( KSPSetType(ksp_, KSPDGMRES) );
  }
  else
    SANS_DEVELOPER_EXCEPTION("UNRECOGNIZED KSP TYPE");


  // Set the preconditioner according to the dictionaries
  setPreconditioner(ksp_, PreconditionerDict_);

  // Set up restart value
  PETSc_STATUS( KSPGMRESSetRestart(ksp_, GMRES_restart_) );
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
void PETScSolver<Matrix_type>::factorize_petsc()
{
  PETSc_STATUS( MatAssemblyBegin(A_petsc_, MAT_FINAL_ASSEMBLY) );
  PETSc_STATUS( MatAssemblyEnd(A_petsc_, MAT_FINAL_ASSEMBLY) );

  // Set the matrix for the linear solver so the preconditioner gets updated
  PETSc_STATUS( KSPSetOperators(ksp_, A_petsc_, A_petsc_) );
  PETSc_STATUS( KSPSetUp(ksp_) );
}


//-----------------------------------------------------------------------------
template< class Matrix_type >
LinearSolveStatus PETScSolver<Matrix_type>::backsolve( const SparseVectorView_type& b, SparseVectorView_type& x ) const
{
  return backsolve(transpose_, b, x);
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
LinearSolveStatus PETScSolver<Matrix_type>::backsolveTranspose( const SparseVectorView_type& b, SparseVectorView_type& x ) const
{
  return backsolve(true, b, x);
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
void
PETScSolver<Matrix_type>::setPreconditioner(KSP ksp, const PyDict& PreconditionerDict)
{
  //Get the preconditioner type from the dictionary
  std::string pc_type = PreconditionerDict.get(params.Preconditioner.Name);

  // Get the preconditioner for the KSP
  PC pc;
  PETSc_STATUS( KSPGetPC(ksp, &pc) );

  std::string pc_side = PreconditionerDict.get(PreconditionerNoneParam::params.PreconditionerSide);
  if (pc_side == PreconditionerNoneParam::params.PreconditionerSide.Left)
  {
    PETSc_STATUS( KSPSetPCSide(ksp, PC_LEFT) );
  }
  else if (pc_side == PreconditionerNoneParam::params.PreconditionerSide.Right)
  {
    PETSc_STATUS( KSPSetPCSide(ksp, PC_RIGHT) );
  }
  else
    SANS_DEVELOPER_EXCEPTION("Unknown preconditioner side.");


  if (pc_type == params.Preconditioner.PETScDefault)
  {
    // do nothing
  }
  else if (pc_type == params.Preconditioner.BlockJacobi)
  {
    PETSc_STATUS( PCSetType(pc, PCBJACOBI) );

    // create the sub-preconditioners
    PETSc_STATUS( KSPSetUp(ksp) );

    // get the sup precondinioner
    PyDict SubPCDict = PreconditionerDict.get(PreconditionerBlockJacobiParam::params.SubPreconditioner);

    PetscInt nlocal,first;
    KSP *subksp;     /* array of local KSP contexts on this processor */

    // Extract the array of KSP contexts for the local blocks
    PETSc_STATUS( PCBJacobiGetSubKSP(pc, &nlocal, &first, &subksp) );

    // Loop over the local blocks, setting various KSP options for each block.
    for (PetscInt i = 0; i < nlocal; i++)
      setPreconditioner(subksp[i], SubPCDict);
  }
  else if (pc_type == params.Preconditioner.ASM)
  {
    PETSc_STATUS( PCSetType(pc, PCASM) );

    PetscInt overlap = PreconditionerDict.get(PreconditionerASMParam::params.Overlap);
    PETSc_STATUS( PCASMSetOverlap(pc, overlap) );

    // create the sub-preconditioners
    PETSc_STATUS( KSPSetUp(ksp) );

    // get the sup precondinioner
    PyDict SubPCDict = PreconditionerDict.get(PreconditionerASMParam::params.SubPreconditioner);

    PetscInt nlocal, first;
    KSP *subksp;     /* array of local KSP contexts on this processor */

    // Extract the array of KSP contexts for the local blocks
    PETSc_STATUS( PCASMGetSubKSP(pc, &nlocal, &first, &subksp) );

    // Loop over the local blocks, setting various KSP options for each block.
    for (PetscInt i = 0; i < nlocal; i++)
    {
      setPreconditioner(subksp[i], SubPCDict);
      KSPSetTolerances(subksp[i], RelativeTolerance_, AbsoluteTolerance_, DivergenceTolerance_, MaxIterations_);

      if (solverType_ == PETScSolverParam::params.KSPSolver.GMRES) {}
      else if (solverType_ == PETScSolverParam::params.KSPSolver.BICGStab)
      {
        PETSc_STATUS( KSPSetType(subksp[i], KSPBCGSL) );
      }
      else if (solverType_ == PETScSolverParam::params.KSPSolver.DGMRES)
      {
        PETSc_STATUS( KSPSetType(subksp[i], KSPDGMRES) );
      }
      else
        SANS_DEVELOPER_EXCEPTION("UNRECOGNIZED KSP TYPE");
    }
  }
  else if (pc_type == params.Preconditioner.Jacobi)
  {
    PETSc_STATUS( PCSetType(pc, PCJACOBI) );
  }
  else if (pc_type == params.Preconditioner.HYPRE)
  {
    PETSc_STATUS( PCSetType(pc, PCHYPRE) );
  }
  else if (pc_type == params.Preconditioner.ILU)
  {
    PETSc_STATUS( PCSetType(pc, PCILU) );

    PetscInt fill_level = PreconditionerDict.get(PreconditionerILUParam::params.FillLevel);
    PETSc_STATUS( PCFactorSetLevels(pc, fill_level) );

    bool nonzero_diagonals = PreconditionerDict.get(PreconditionerILUParam::params.NonzeroDiagonal);
    if (nonzero_diagonals)
      PETSc_STATUS( PCFactorReorderForNonzeroDiagonal(pc, 1000*std::numeric_limits<Real>::epsilon()) );

    std::string ordering = PreconditionerDict.get(PreconditionerILUParam::params.Ordering);
    if (ordering == PreconditionerILUParam::params.Ordering.Natural)
    {
      PETSc_STATUS( PCFactorSetMatOrderingType(pc, MATORDERINGNATURAL) );
    }
    else if (ordering == PreconditionerILUParam::params.Ordering.NestedDissection)
    {
      PETSc_STATUS( PCFactorSetMatOrderingType(pc, MATORDERINGND) );
    }
    else if (ordering == PreconditionerILUParam::params.Ordering.ReverseCuthillMcKee)
    {
      PETSc_STATUS( PCFactorSetMatOrderingType(pc, MATORDERINGRCM) );
    }
    else if (ordering == PreconditionerILUParam::params.Ordering.QDM)
    {
      PETSc_STATUS( PCFactorSetMatOrderingType(pc, MATORDERINGQMD) );
    }
    else if (ordering == PreconditionerILUParam::params.Ordering.MDF)
    {
      PETSc_STATUS( MatOrderingRegister("mdf", MatGetOrdering_MDF) );
      PETSc_STATUS( PCFactorSetMatOrderingType(pc, "mdf") );

      MDFOuterBlockSize_ = PreconditionerDict.get(PreconditionerILUParam::params.MDFOuterBlockSize);

      //Hack to send outerblocksize to MDF routine (which is static)
      PETSc_STATUS( PetscContainerSetPointer(mdfcontainer_, &MDFOuterBlockSize_ ) );

      //Compose the preconditioner matrix with a pointer to the MDFouterblocksize, so that it can be queried from the static ordering function
      Mat mat;
      PETSc_STATUS( PCGetOperators(pc, &mat, NULL) );
      PETSc_STATUS( PetscObjectCompose( (PetscObject) mat, "MDFouterblocksize", (PetscObject) mdfcontainer_) );
      PETSc_STATUS( PetscObjectCompose( (PetscObject) A_petsc_, "MDFouterblocksize", (PetscObject) mdfcontainer_) );
    }
    else
      SANS_DEVELOPER_EXCEPTION("Unknown ordering type for ILU preconditioner.");
  }
  else if (pc_type == params.Preconditioner.DirectLU)
  {
    PETSc_STATUS( PCSetType(pc, PCLU) );
#ifdef PETSC_HAVE_MKL_PARDISO
    //PETSc_STATUS( PCFactorSetMatSolverType(pc, MATSOLVERMKL_PARDISO) );
#endif
  }
  else if (pc_type == params.Preconditioner.None)
  {
    PETSc_STATUS( PCSetType(pc, PCNONE) );
  }
  else
    SANS_DEVELOPER_EXCEPTION("Unknown preconditioner type.");
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
void
PETScSolver<Matrix_type>::view(PetscViewer viewer) const
{
  if (viewer == nullptr)
    viewer =  PETSC_VIEWER_STDOUT_(*continuousmap_[0].comm);

  PETSc_STATUS( MatView(A_petsc_, viewer) );
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
void
PETScSolver<Matrix_type>::getOrderingIndex(MatOrderingType type, IS *rperm, IS *cperm)
{
  MatType mattype;
  PETSc_STATUS( MatGetType(A_petsc_, &mattype) );

  // Get the preconditioner for the KSP
  PC pc;
  PETSc_STATUS( KSPGetPC(ksp_, &pc) );

  if ( strcmp(mattype, MATSEQAIJ) == 0 || strcmp(mattype, MATMPIAIJ) == 0 || strcmp(mattype, MATSEQBAIJ) == 0 )
  {
    Mat mat;
    PETSc_STATUS( PCGetOperators(pc, &mat, NULL) );

    //These matrix types work fine with MatGetOrdering
    PETSc_STATUS( MatGetOrdering(mat, type, rperm, cperm) );
    return;
  }
  else if ( strcmp(mattype, MATMPIBAIJ) == 0 )
  {
    //Calling MatGetOrdering directly using A_petsc_ only returns a natural ordering,
    //so we need to extract the local sub-matrix and call individually from each processor

    // Get the local submatrix on this processor
    PetscInt nMat = 1; //default is one local matrix per processor
    Mat *matLocal;
    PETSc_STATUS( PCASMGetLocalSubmatrices(pc, &nMat, &matLocal) );
    PETSc_STATUS( MatGetOrdering(matLocal[0], type, rperm, cperm) );
  }
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
void
PETScSolver<Matrix_type>::getDebugInfo() const
{
  if (verbose_ && continuousmap_[0].comm->rank() == 0)
  {
    KSPConvergedReason reason;
    PetscInt iter = -1;
    PetscReal rnorm = -1;
    PETSc_STATUS( KSPGetConvergedReason(ksp_, &reason) );
    PETSc_STATUS( KSPGetIterationNumber(ksp_, &iter) );
    PETSc_STATUS( KSPGetResidualNorm(ksp_, &rnorm) );
    std::cout << "PETSc " << iter << ", " << std::scientific << std::setprecision(5) << rnorm << std::endl;
    std::cout << "PETSc result: " << getReasonDescription(reason) << std::endl;
  }

  if (printMatrixInfo_)
  {
    MatInfo info;
    PetscLogDouble nz_allocated = 0, nz_used = 0, nz_unneeded = 0;
    PetscLogDouble memory = 0;
    PETSc_STATUS( MatGetInfo(A_petsc_, MAT_LOCAL, &info) );

#ifdef SANS_MPI
    nz_allocated = boost::mpi::all_reduce( *continuousmap_[0].comm, info.nz_allocated, std::plus<Real>() );
    nz_used = boost::mpi::all_reduce( *continuousmap_[0].comm, info.nz_used, std::plus<Real>() );
    nz_unneeded = boost::mpi::all_reduce( *continuousmap_[0].comm, info.nz_unneeded, std::plus<Real>() );
    memory = boost::mpi::all_reduce( *continuousmap_[0].comm, info.memory, std::plus<Real>() );
#else
    nz_allocated = info.nz_allocated;
    nz_used = info.nz_used;
    nz_unneeded = info.nz_unneeded;
    memory = info.memory;
#endif

    if (continuousmap_[0].comm->rank() == 0)
    {
      std::cout << "PETSc matrix info:" << std::endl;
      std::cout << "  Nonzeros allocated: " << nz_allocated << std::endl;
      std::cout << "  Nonzeros used: " << nz_used << std::endl;
      std::cout << "  Nonzeros unneeded: " << nz_unneeded << std::endl;
      std::cout << "  Memory used: " << memory/1.0e6 << " MB" << std::endl;
    }
  }

  if (computeSingularValues_ && continuousmap_[0].comm->rank() == 0)
  {
    PetscReal sMax, sMin;
    PETSc_STATUS( KSPComputeExtremeSingularValues(ksp_, &sMax, &sMin) );
    std::cout << "PETSc sigmaMax: " << sMax << std::endl;
    std::cout << "PETSc sigmaMin: " << sMin << std::endl;
    std::cout << "PETSc sigmaMax/sigmaMin: " << sMax/sMin << std::endl;

#if 0 //If we need to compute multiple extreme eigenvalues
    int n = 50;
    std::vector<PetscReal> eigRe(n), eigIm(n);
    PetscInt nEig = 0;
    PETSc_STATUS( KSPComputeEigenvalues(ksp_, n, eigRe.data(), eigIm.data(), &nEig) );

    for (int i = 0; i < nEig; i++)
    {
      PetscReal eigAbs = sqrt(eigRe[i]*eigRe[i] + eigIm[i]*eigIm[i]);
      std::cout << "Eig" << i << ": " << eigRe[i] << " + " << eigIm[i] << "i : " << eigAbs << std::endl;
    }
#endif
  }

#if 0
  Vec res;
  PETSc_STATUS( KSPBuildResidual(ksp_, NULL, NULL, &res) );

  PetscInt m;
  PETSc_STATUS( VecGetSize(res, &m) );

  std::string filename = filename_base_ + "PETSc_resvec.dat";
  std::fstream fresidual( filename, std::fstream::out );
  fresidual << std::scientific << std::setprecision(10);

  PetscReal val;
  for (int i = 0; i < m; i++)
  {
    PETSc_STATUS( VecGetValues(res, 1, &i, &val) );
    fresidual << val << std::endl;
  }
  fresidual.close();
#endif
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
std::string PETScSolver<Matrix_type>::getReasonDescription(KSPConvergedReason& reason) const
{
  switch (reason)
  {
    case KSP_CONVERGED_RTOL_NORMAL:
      return "Converged - Relative tolerance normal";
      break;

    case KSP_CONVERGED_ATOL_NORMAL:
      return "Converged - Absolute tolerance normal";
      break;

    case KSP_CONVERGED_RTOL:
      return "Converged - Residual 2-norm decreased by a factor of rtol";
      break;

    case KSP_CONVERGED_ATOL:
      return "Converged - Residual 2-norm less than abstol";
      break;

    case KSP_DIVERGED_NULL:
      return "Diverged - Null";
      break;

    case KSP_DIVERGED_ITS:
      return "Diverged - Maximum iteration count reached";
      break;

    case KSP_DIVERGED_DTOL:
      return "Diverged - Residual norm increased by a factor of divtol";
      break;

    case KSP_DIVERGED_NANORINF:
      return "Diverged - Residual norm became Not-a-Number or Inf likely due to 0/0";
      break;

    case KSP_DIVERGED_BREAKDOWN:
      return "Diverged - Generic breakdown in method";
      break;

#if PETSC_VERSION_LT(3,11,0)
    case KSP_DIVERGED_PCSETUP_FAILED:
#else
    case KSP_DIVERGED_PC_FAILED:
#endif
      return "Diverged - Preconditioner failed";
      break;

    case KSP_CONVERGED_ITERATING:
      return "KSP Solve is still running!";
      break;

    default:
      return "Code " + std::to_string(reason) + " - Unknown description";
      break;
  }
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
PetscErrorCode
PETScSolver<Matrix_type>::residualMonitor(KSP ksp, PetscInt it, PetscReal rnorm, void *mctx)
{
  ResidualMonitor& rsdMon = *static_cast<ResidualMonitor*>(mctx);

  if ( rsdMon.verbose && (rsdMon.comm_rank == 0) && (it % (rsdMon.MaxIterations/5) == 0) )
    std::cout << "PETSc " << it << ", " << std::scientific << std::setprecision(5) << rnorm << std::endl;

  if ( rsdMon.dumpResidualHistory && (rsdMon.comm_rank == 0) )
    rsdMon.fhist << it << ", " << std::scientific << std::setprecision(16) << rnorm << std::endl;

  return 0;
}

//-----------------------------------------------------------------------------
template< class Matrix_type >
PETScSolver<Matrix_type>::~PETScSolver()
{
  delete this->A_;
  delete pMs_;
  PETSc_STATUS( KSPDestroy(&ksp_) );
  PETSc_STATUS( MatDestroy(&A_petsc_) );
  PETSc_STATUS( PetscContainerDestroy(const_cast<PetscContainer*>(&mdfcontainer_)) );
}

} //namespace SLA
} //namespace numpack 
