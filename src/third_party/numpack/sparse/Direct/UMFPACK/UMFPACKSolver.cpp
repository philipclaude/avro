// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "UMFPACKSolver.h"
#include "tools/stringify.h"

#include <umfpack.h>

namespace numpack 
{
namespace SLA
{

//=============================================================================
UMFPACKException::UMFPACKException( const int status )
{
  errString += "UMFPACK Error\n";

  switch (status)
  {
  case UMFPACK_ERROR_invalid_matrix:
    errString += "UMFPACK: Invalid matrix\n";
    break;
  case UMFPACK_ERROR_out_of_memory:
    errString += "UMFPACK: Out of Memory\n";
    break;
  case UMFPACK_ERROR_internal_error:
    errString += "UMFPACK: Internal Error. You found a bug in UMFPACK!.\nPlease contact the author (DrTimothyAldenDavis@gmail.com).\n";
    break;
  case UMFPACK_WARNING_singular_matrix:
    errString += "UMFPACK: Singular Matrix!\n";
    break;
  case UMFPACK_ERROR_invalid_Symbolic_object:
    errString += "UMFPACK: Invalid Symbolic Object\n";
    break;
  case UMFPACK_ERROR_different_pattern:
    errString += "UMFPACK: The pattern has changed. This is a developer error.\n";
    break;
  }

}

//=============================================================================
UMFPACKException::UMFPACKException( const int status, const double* info )
{
  errString += "UMFPACK Error'\n";

  double Unit = info[UMFPACK_SIZE_OF_UNIT];

  switch (status)
  {
  case UMFPACK_ERROR_invalid_matrix:
    errString += "UMFPACK: Invalid matrix\n";
    break;
  case UMFPACK_ERROR_out_of_memory:
    errString += "UMFPACK: Out of Memory!\n";
    errString += "Symbolic memory usage: " + stringify(info[UMFPACK_SYMBOLIC_PEAK_MEMORY]*Unit) + " bytes\n";
    errString += "Numeric  memory usage: " + stringify(info[UMFPACK_NUMERIC_SIZE]*Unit) + " bytes\n";
    errString += "Total    memory usage: " + stringify(info[UMFPACK_PEAK_MEMORY]*Unit) + " bytes\n";
    break;
  case UMFPACK_ERROR_internal_error:
    errString += "UMFPACK: Internal Error. You found a bug in UMFPACK!.\nPlease contact the author (DrTimothyAldenDavis@gmail.com).\n";
    break;
  case UMFPACK_WARNING_singular_matrix:
    errString += "UMFPACK: Singular Matrix!\n";
    break;
  case UMFPACK_ERROR_invalid_Symbolic_object:
    errString += "UMFPACK: Invalid Symbolic Object\n";
    break;
  case UMFPACK_ERROR_different_pattern:
    errString += "UMFPACK: The pattern has changed. This is a developer error.\n";
    break;
  }
}

// cppcheck-suppress passedByValue
void UMFPACKParam::checkInputs(PyDict d)
{
  std::vector<const ParameterBase*> allParams;
  allParams.push_back(d.checkInputs(params.Timing));
  d.checkUnknownInputs(allParams);
}
UMFPACKParam UMFPACKParam::params;


} //namespace SLA
} //namespace numpack 
