// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef PETSCSOLVER_H
#define PETSCSOLVER_H

#include "Python/PyDict.h"
#include "Python/Parameter.h"

#include "tools/SANSException.h"

#include "numpack/sparse/LinearSolverBase.h"
#include "numpack/sparse/sparse_Inverse.h"

#include "numpack/AlgebraicEquationSetBase.h"
#include "numpack/GlobalContinuousMap.h"

#include <memory>
#include <fstream>

#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>
#include <petscpc.h>

namespace numpack 
{
namespace SLA
{

//=============================================================================
struct PETScSolverException : public SANSException
{
  explicit PETScSolverException(const PetscErrorCode status);

  virtual ~PETScSolverException() throw() {}
};

#define PETSc_STATUS( API_call ) \
  {  \
    PetscErrorCode _eg_status_from_api_call_ = API_call; \
    if ( _eg_status_from_api_call_ != 0 ) \
      BOOST_THROW_EXCEPTION( numpack::SLA::PETScSolverException(_eg_status_from_api_call_) ); \
  }

//=============================================================================
//Forward declare
template< class Matrix_type >
class PETScSolver;

//=============================================================================

struct PreconditionerCommonParam : noncopyable
{
  struct PreconditionerSideOptions
  {
    typedef std::string ExtractType;
    const std::string Left = "Left";
    const std::string Right = "Right";

    const std::vector<std::string> options{Left, Right};
  };
  const ParameterOption<PreconditionerSideOptions> PreconditionerSide =
      ParameterOption<PreconditionerSideOptions>("PreconditionerSide", "Right", "Preconditioning side");

  struct OrderingOptions
  {
    typedef std::string ExtractType;
    const std::string Natural = "Natural";
    const std::string NestedDissection = "NestedDissection";
    const std::string ReverseCuthillMcKee = "ReverseCuthillMcKee";
    const std::string QDM = "QDM";
    const std::string MDF = "MDF";

    const std::vector<std::string> options{Natural, NestedDissection, ReverseCuthillMcKee, QDM, MDF};
  };
  const ParameterOption<OrderingOptions> Ordering =
      ParameterOption<OrderingOptions>("Ordering", "Natural", "Matrix ordering type");

  // This is a hack to make the MDF ordering operate on elemental blocks instead of PDE blocks
  // MDFOuterBlockSize should be set to the number of DOFs in each element (e.g. 3 for DG P1 triangles)
  const ParameterNumeric<int> MDFOuterBlockSize{"MDFOuterBlockSize", 1, 1, NO_LIMIT, "Size of outer block (MatrixD) for MDF ordering"};
};

struct PreconditionerILUParam : public PreconditionerCommonParam
{
  const ParameterNumeric<int> FillLevel{"FillLevel", 0, 0, 10, "Fill-level - ILU(0), ILU(1) etc"};
  const ParameterBool NonzeroDiagonal{"NonzeroDiagonal", false, "Re-order for non-zero diagonals?"};

  static void checkInputs(PyDict d);
  static PreconditionerILUParam params;
};

struct PreconditionerJacobiParam : public PreconditionerCommonParam
{
  static void checkInputs(PyDict d);
  static PreconditionerJacobiParam params;
};

struct PreconditionerHYPREParam : public PreconditionerCommonParam
{
  static void checkInputs(PyDict d);
  static PreconditionerHYPREParam params;
};

struct PreconditionerDirectLUParam : public PreconditionerCommonParam
{
  static void checkInputs(PyDict d);
  static PreconditionerDirectLUParam params;
};

struct PreconditionerNoneParam : public PreconditionerCommonParam
{
  static void checkInputs(PyDict d);
  static PreconditionerNoneParam params;
};

struct SubPreconditionerParam : public PreconditionerCommonParam
{
  struct PreconditionerOptions
  {
    typedef DictKeyPair ExtractType;
    const ParameterString Name{"Name", "PETScDefault", "Sub-Preconditioner name" };
    const ParameterString& key = Name;

    const DictOption PETScDefault{"PETScDefault", PreconditionerNoneParam::checkInputs};
    const DictOption Jacobi{"Jacobi", PreconditionerJacobiParam::checkInputs};
    const DictOption ILU{"ILU", PreconditionerILUParam::checkInputs};
    const DictOption DirectLU{"DirectLU", PreconditionerDirectLUParam::checkInputs};
    const DictOption None{"None", PreconditionerNoneParam::checkInputs};

    const std::vector<DictOption> options{PETScDefault, Jacobi, ILU, DirectLU, None};
  };
  const ParameterOption<PreconditionerOptions> SubPreconditioner{"SubPreconditioner", EMPTY_DICT, "Sub Preconditioner"};

  struct PreconditionerSideOptions
  {
    typedef std::string ExtractType;
    const std::string Left = "Left";
    const std::string Right = "Right";

    const std::vector<std::string> options{Left, Right};
  };
  const ParameterOption<PreconditionerSideOptions> PreconditionerSide =
      ParameterOption<PreconditionerSideOptions>("PreconditionerSide", "Left", "Preconditioning side");
};

struct PreconditionerBlockJacobiParam : public SubPreconditionerParam
{
  static void checkInputs(PyDict d);
  static PreconditionerBlockJacobiParam params;
};

struct PreconditionerASMParam : public SubPreconditionerParam
{
  const ParameterNumeric<int> Overlap{"Overlap", 1, 0, NO_LIMIT, "ASM Overlap on subdomains"};

  static void checkInputs(PyDict d);
  static PreconditionerASMParam params;
};

struct PETScSolverParam : noncopyable
{

  template<class Matrix_type>
  static std::shared_ptr< LinearSolverBase<Matrix_type> >
  newSolver(PyDict& SolverParam, AlgebraicEquationSetBase<Matrix_type>& f, LinearSystemSolve solve = RegularSolve)
  {
    typedef std::shared_ptr< LinearSolverBase<Matrix_type> > Solver_ptr;

    return Solver_ptr( new PETScSolver<Matrix_type>( SolverParam, f, solve ) );
  }

  const ParameterNumeric<Real> RelativeTolerance{"RelativeTolerance", 1e-05, 0, NO_LIMIT, "Relative Convergence Tolerance "};
  const ParameterNumeric<Real> AbsoluteTolerance{"AbsoluteTolerance", 1e-50, 0, NO_LIMIT, "Absolute Convergence Tolerance"};
  const ParameterNumeric<Real> DivergenceTolerance{"DivergenceTolerance", 1e7, 1, NO_LIMIT, "Divergence Tolerance"};
  const ParameterNumeric<int > MaxIterations{"MaxIterations", 10000, 0, NO_LIMIT, "Maximum Number of Linear Solver Iterations"};
  const ParameterNumeric<int > GMRES_Restart{"GMRES_Restart", 30, 1, NO_LIMIT, "Number of iterations after which GMRES restarts"};

  struct PreconditionerOptions
  {
    typedef DictKeyPair ExtractType;
    const ParameterString Name{"Name", "PETScDefault", "Preconditioner name" };
    const ParameterString& key = Name;

    const DictOption PETScDefault{"PETScDefault", PreconditionerNoneParam::checkInputs};
    const DictOption BlockJacobi{"BlockJacobi", PreconditionerBlockJacobiParam::checkInputs};
    const DictOption Jacobi{"Jacobi", PreconditionerILUParam::checkInputs};
    const DictOption ILU{"ILU", PreconditionerILUParam::checkInputs};
    const DictOption ASM{"ASM", PreconditionerASMParam::checkInputs};
    const DictOption HYPRE{"HYPRE", PreconditionerHYPREParam::checkInputs};
    const DictOption DirectLU{"DirectLU", PreconditionerDirectLUParam::checkInputs};
    const DictOption None{"None", PreconditionerNoneParam::checkInputs};

    const std::vector<DictOption> options{PETScDefault, BlockJacobi, Jacobi, ILU, ASM, HYPRE, DirectLU, None};
  };


  struct KSPSolverOptions
  {
    typedef std::string ExtractType;
    const std::string GMRES = "GMRES";
    const std::string DGMRES = "DGMRES"; //deflated gmres
    const std::string BICGStab = "BICGStab";

    const std::vector<std::string> options{GMRES, DGMRES, BICGStab};
  };
  const ParameterOption<KSPSolverOptions> KSPSolver =
      ParameterOption<KSPSolverOptions>("KSPSolver", "GMRES", "KSP Solver Type");

  const ParameterOption<PreconditionerOptions> Preconditioner{"Preconditioner", EMPTY_DICT, "Preconditioner"};

  const ParameterBool Verbose{"Verbose", false, "Verbose output?"};
  const ParameterBool Timing{"Timing", false, "Time Components of PETSc Solve"};
  const ParameterBool Memory{"Memory", false, "Report Memory Usage"};
  const ParameterBool computeSingularValues{"computeSingularValues", false, "Compute extreme singular values?"};
  const ParameterBool printMatrixInfo{"printMatrixInfo", false, "Print matrix info (nnz)?"};
  const ParameterFileName ResidualHistoryFile{"ResidualHistoryFile", std::ios_base::out, "", "Dumps residual history to give file name"};
  const ParameterString FilenameBase{"FilenameBase", "tmp/", "Default filepath prefix for debug files"};

  const ParameterString PETScOptions{"PETScOptions", "", "String of PETSc command line options"};
  const ParameterString PETScOptionsFile{"PETScOptionsFile", "", "File name of PETSc command line options"};

  const ParameterString PETScMatPrefix{"PETScMatPrefix", "", "Prefix for PETSc Mat command line options"};
  const ParameterString PETScKSPPrefix{"PETScKSPPrefix", "", "Prefix for PETSc KSP command line options"};

  static void checkInputs(PyDict d);
  static PETScSolverParam params;
};


//=============================================================================
template< class Matrix_type >
class PETScSolver : public LinearSolverBase< Matrix_type >
{
public:
  typedef LinearSolverBase< Matrix_type > Base_type;

  typedef typename Base_type::SparseVector_type SparseVector_type;
  typedef typename Base_type::SparseVectorView_type SparseVectorView_type;
  typedef typename NonZeroPatternType<Matrix_type>::type SystemNonZeroPattern;

  PETScSolverParam& params;

//-----------------------------------------------------------------------------
  PETScSolver(PyDict d, AlgebraicEquationSetBase<Matrix_type>& f, LinearSystemSolve solve = RegularSolve);
  virtual ~PETScSolver();

//-----------------------------------------------------------------------------
  virtual void factorize() override;

  //-----------------------------------------------------------------------------
  void factorize( SparseVectorView_type& bcondensed, bool transpose); //for static condensation

//-----------------------------------------------------------------------------
  virtual LinearSolveStatus backsolve( const SparseVectorView_type& b, SparseVectorView_type& x ) const override;

//-----------------------------------------------------------------------------
  LinearSolveStatus backsolveTranspose( const SparseVectorView_type& b, SparseVectorView_type& x ) const;

//-----------------------------------------------------------------------------
  virtual LinearSolveStatus solve(SparseVectorView_type& b, SparseVectorView_type& x) override;

//-----------------------------------------------------------------------------
  const KSP& getKSP() const { return ksp_; }
  const Mat& getMat() const { return A_petsc_; };

  void view(PetscViewer viewer = nullptr) const;

  void getOrderingIndex(MatOrderingType type, IS *rperm, IS *cperm);

protected:
//-----------------------------------------------------------------------------
  void init();
  void init_petsc();

  void factorizeMatrix();
  void factorize_petsc();

  LinearSolveStatus backsolve( const bool transpose, const SparseVectorView_type& b, SparseVectorView_type& x ) const;

//-----------------------------------------------------------------------------
  void setPreconditioner(KSP ksp, const PyDict& PreconditionerDict);

  static PetscErrorCode MatGetOrdering_MDF(Mat mat, MatOrderingType type, IS *irow, IS *icol);

  void getDebugInfo() const;

  std::string getReasonDescription(KSPConvergedReason& reason) const;

  static PetscErrorCode residualMonitor(KSP ksp, PetscInt it, PetscReal rnorm, void *mctx);

protected:

  struct ResidualMonitor
  {
    ResidualMonitor(std::string filename, bool verbose, PetscInt MaxIterations, int comm_rank);

    void setCommRank(const int comm_rank);

    const bool dumpResidualHistory;
    std::ofstream fhist;
    const bool verbose;
    const PetscInt MaxIterations;
    const int comm_rank;
  };

  using Base_type::A_;
  using Base_type::transpose_;
  AlgebraicEquationSetBase<Matrix_type>& f_;
  std::vector<GlobalContinuousMap> continuousmap_;
  ScalarMatrix_CRS<PetscInt, PetscScalar> *pMs_;
  Mat A_petsc_;          // linear system matrix
  KSP ksp_;              // linear solver context
  // TBD
  PetscReal RelativeTolerance_;
  PetscReal AbsoluteTolerance_;
  PetscReal DivergenceTolerance_;
  PetscInt MaxIterations_;
  PetscInt GMRES_restart_;
  PyDict PreconditionerDict_;
  bool verbose_;
  bool timing_;
  bool memory_;
  bool computeSingularValues_;
  bool printMatrixInfo_;
  std::string matprefix_;
  std::string kspprefix_;
  ResidualMonitor rsdMon_;
  PetscInt MDFOuterBlockSize_;
  PetscContainer mdfcontainer_;
  std::string solverType_;
  bool staticCondensed_;
};

} //namespace SLA
} //namespace numpack 

#endif //PETSCSOLVER_H
