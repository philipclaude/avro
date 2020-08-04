// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef ALGEBRAICEQUATIONSETBASE_H
#define ALGEBRAICEQUATIONSETBASE_H

#include "tools/SANSnumerics.h"
#include "tools/SANSException.h"

#include <vector>
#include <string>
#include <memory> // std::shared_ptr
#include <iostream> // std:cout

#include "VectorType.h"
#include "MatrixViewType.h"
#include "NonZeroPatternType.h"
#include "VectorSizeType.h"
#include "MatrixSizeType.h"

//#include "MPI/communicator_fwd.h"
#include "LinesearchDataType.h"
#include "GlobalContinuousMap.h"

namespace numpack
{

//----------------------------------------------------------------------------//
// An abstract base class for residual operations
//
//
// template parameters:
//

template<class SystemMatrix_>
class AlgebraicEquationSetBase
{
public:
  typedef SystemMatrix_ SystemMatrix;
  typedef typename VectorType<SystemMatrix>::type SystemVector;
  typedef typename NonZeroPatternType<SystemMatrix>::type SystemNonZeroPattern;

  typedef typename MatrixViewType<SystemMatrix>::type SystemMatrixView;
  typedef typename VectorType<SystemMatrix>::Viewtype SystemVectorView;
  typedef typename NonZeroPatternType<SystemMatrix>::Viewtype SystemNonZeroPatternView;

  typedef typename MatrixSizeType<SystemMatrix>::type MatrixSizeClass;
  typedef typename VectorSizeType<SystemVector>::type VectorSizeClass;

  typedef typename LinesearchDataType<SystemVector>::type LinesearchData;

  virtual ~AlgebraicEquationSetBase() {}

  //Computes the residual
  virtual void residual(                           SystemVectorView& rsd) = 0;
          void residual(const SystemVectorView& q, SystemVectorView& rsd) { setSolutionField(q); residual(rsd); }

  //Computes the jacobian or fills the non-zero pattern of a jacobian
  virtual void jacobian(                           SystemMatrixView& mtx       ) = 0;
  virtual void jacobian(                           SystemNonZeroPatternView& nz) = 0;
          void jacobian(const SystemVectorView& q, SystemMatrixView& mtx       ) { setSolutionField(q); jacobian(mtx); }
          void jacobian(const SystemVectorView& q, SystemNonZeroPatternView& nz) { setSolutionField(q); jacobian(nz);  }
  virtual void jacobianTranspose(                           SystemMatrixView& mtx       ) = 0;
  virtual void jacobianTranspose(                           SystemNonZeroPatternView& nz) = 0;
          void jacobianTranspose(const SystemVectorView& q, SystemMatrixView& mtx       ) { setSolutionField(q); jacobianTranspose(mtx); }
          void jacobianTranspose(const SystemVectorView& q, SystemNonZeroPatternView& nz) { setSolutionField(q); jacobianTranspose(nz);  }

  virtual void jacobian(SystemVectorView& rsd, SystemNonZeroPatternView& nz, bool transpose )
  {
    SANS_DEVELOPER_EXCEPTION("SHOULDN'T BE USING THIS NZ CALL UNLESS IMPLEMENTED FOR STATIC CONDENSATION");
  }

  virtual void jacobian(SystemVectorView& rsd, SystemMatrixView& mtx, bool transpose )
  {
    SANS_DEVELOPER_EXCEPTION("SHOULDN'T BE USING THIS JACOBIAN CALL UNLESS IMPLEMENTED FOR STATIC CONDENSATION");
  }

  //Compute Residual Norm
  virtual std::vector<std::vector<Real>> residualNorm(const SystemVectorView& rsd) const = 0;

  //Convergence check of the residual
  virtual bool convergedResidual(const std::vector<std::vector<Real>>& rsdNorm) const = 0;
  virtual bool convergedResidual(const std::vector<std::vector<Real>>& rsdNorm, int iEq, int iMon) const = 0;

  //Check whether residual decreased, for purposes of linesearch
  virtual bool decreasedResidual(const std::vector<std::vector<Real>>& rsdNormOld,
                                 const std::vector<std::vector<Real>>& rsdNormNew) const = 0;

  //prints out a residual that could not be decreased and the convergence tolerances
  virtual void printDecreaseResidualFailure(const std::vector<std::vector<Real>>& rsdNorm, std::ostream& os = std::cout) const = 0;

  //Translates the system vector into a solution field
  virtual void setSolutionField(const SystemVectorView& q) = 0;

  //Translates the solution field into a system vector
  virtual void fillSystemVector(SystemVectorView& q) const = 0;

  //Returns the vector and matrix sizes needed for the linear algebra system
  virtual VectorSizeClass vectorEqSize() const = 0;    // vector for equations (rows in matrixSize)
  virtual VectorSizeClass vectorStateSize() const = 0; // vector for state DOFs (columns in matrixSize)
  virtual MatrixSizeClass matrixSize() const = 0;

  virtual bool isStaticCondensed() const { return isStaticCondensed_;}

  virtual void completeUpdate(const SystemVectorView& rsd, SystemVectorView& xcondensed, SystemVectorView& x) const
  {
    SANS_DEVELOPER_EXCEPTION("SHOULDN'T BE USING COMPLETEUPDATE UNLESS IMPLEMENTED FOR STATIC CONDENSATION");
  }

  //Gives the indices to the equations and solution variable blocks in the system
  //Everything is set to -1 in the base class, derived classes should override appropriately
  virtual int indexPDE() const { return -1; }
  virtual int indexBC()  const { return -1; }
  virtual int indexAUX() const { return -1; }
  virtual int indexINT() const { return -1; }
  virtual int indexQ()   const { return -1; }

  virtual int indexPDESC()   const { return -1; }
  virtual int indexQSC()   const { return -1; }

  // update fraction needed for physically valid state
  virtual Real updateFraction(const SystemVectorView& q, const SystemVectorView& dq, const Real maxChangeFraction) const
  {
    SANS_DEVELOPER_EXCEPTION("Not implemented...");
    return 1;
  }

  //Checks to see if proposed solution is physical
  virtual bool isValidStateSystemVector(SystemVectorView& q) = 0;

  //Returns the size of the residual norm outer vector
  virtual int nResidNorm() const = 0;

  // Converts processor local indexing to a processor continuous indexing
  virtual std::vector<GlobalContinuousMap> continuousGlobalMap() const { return {}; } // TODO: This should be abstract

  // MPI communicator for this algebraic equation set
  //virtual std::shared_ptr<mpi::communicator> comm() const = 0;

  // Syncronize gost/zombie DOFs in all fields that need it
  virtual void syncDOFs_MPI() = 0;

  //Update the localized linesearch parameter info (used for debugging purposes)
  virtual bool updateLinesearchDebugInfo(const Real& s, const SystemVectorView& rsd,
                                         LinesearchData& pResData,
                                         LinesearchData& pStepData) const = 0;

  virtual void dumpLinesearchDebugInfo(const std::string& filenamebase, const int& nonlinear_iter,
                                       const LinesearchData& pStepData) const = 0;

  //Experimental debugging utility
  virtual void dumpSolution(const std::string& filename) const { SANS_ASSERT(false); }

  // Are we having machine precision issues
  virtual bool atMachinePrecision(const SystemVectorView& q, const std::vector<std::vector<Real>>& R0norm) = 0;

protected:
  bool isStaticCondensed_ = false;

};


} //namespace numpack

#endif //ALGEBRAICEQUATIONSETBASE_H
