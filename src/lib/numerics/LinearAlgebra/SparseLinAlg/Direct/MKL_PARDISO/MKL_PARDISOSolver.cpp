// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "MKL_PARDISOSolver.h"
#include "tools/stringify.h"

namespace SANS
{
namespace SLA
{

//=============================================================================
MKL_PARDISOException::MKL_PARDISOException( const int status )
{
  errString += "MKL_PARDISO Error\n";

  switch (status)
  {
  case 42:
    errString += "SANS must be linked with Intel(R) MKL libraries to use the MKL PARDISO solver";
    break;
  case -1:
    errString += "input inconsistent";
    break;
  case -2:
    errString += "not enough memory";
    break;
  case -3:
    errString += "reordering problem";
    break;
  case -4:
    errString += "zero pivot, numerical factorization or iterative refinement problem";
    break;
  case -5:
    errString += "unclassified (internal) error";
    break;
  case -6:
    errString += "reordering failed (matrix types 11 and 13 only)";
    break;
  case -7:
    errString += "diagonal matrix is singular";
    break;
  case -8:
    errString += "32-bit integer overflow problem";
    break;
  case -9:
    errString += "not enough memory for OOC";
    break;
  case -10:
    errString += "error opening OOC files";
    break;
  case -11:
    errString += "read/write error with OOC files";
    break;
  case -12:
    errString += "(pardiso_64 only) pardiso_64 called from 32-bit library";
    break;
  default:
    errString += "Unkown error: status = " + stringify(status);
    break;
  }

}

// cppcheck-suppress passedByValue
void MKL_PARDISOParam::checkInputs(PyDict d)
{
  std::vector<const ParameterBase*> allParams;
  allParams.push_back(d.checkInputs(params.Timing));
  d.checkUnknownInputs(allParams);
}
MKL_PARDISOParam MKL_PARDISOParam::params;


} //namespace SLA
} //namespace SANS
