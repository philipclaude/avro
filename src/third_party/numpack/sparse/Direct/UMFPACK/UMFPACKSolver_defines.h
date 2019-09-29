// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef UMFPACKSOLVER_DEFINES_H
#define UMFPACKSOLVER_DEFINES_H

// UMFPACK interface defines

#include <umfpack.h>

namespace numpack 
{
namespace SLA
{

#define USE_SANS_UMFPACK_LONG

#ifdef USE_SANS_UMFPACK_LONG            // UMFPACK 64-bit interface

  #ifndef SuiteSparse_long
  #define SuiteSparse_long UF_long // Use the depricated UF_long if SuiteSparse_long is not defined
  #endif

  #define SANS_UMFPACK_INT              SuiteSparse_long
  #define SANS_UMFPACK_DEFAULTS         umfpack_dl_defaults
  #define SANS_UMFPACK_FREE_NUMERIC     umfpack_dl_free_numeric
  #define SANS_UMFPACK_FREE_SYMBOLIC    umfpack_dl_free_symbolic
  #define SANS_UMFPACK_NUMERIC          umfpack_dl_numeric
  #define SANS_UMFPACK_REPORT_INFO      umfpack_dl_report_info
  #define SANS_UMFPACK_SOLVE            umfpack_dl_solve
  #define SANS_UMFPACK_SYMBOLIC         umfpack_dl_symbolic
  #define SANS_UMFPACK_TRANSPOSE        umfpack_dl_transpose

#else                                   // UMFPACK 32-bit interface
  #define SANS_UMFPACK_INT              int
  #define SANS_UMFPACK_DEFAULTS         umfpack_di_defaults
  #define SANS_UMFPACK_FREE_NUMERIC     umfpack_di_free_numeric
  #define SANS_UMFPACK_FREE_SYMBOLIC    umfpack_di_free_symbolic
  #define SANS_UMFPACK_NUMERIC          umfpack_di_numeric
  #define SANS_UMFPACK_REPORT_INFO      umfpack_di_report_info
  #define SANS_UMFPACK_SOLVE            umfpack_di_solve
  #define SANS_UMFPACK_SYMBOLIC         umfpack_di_symbolic
  #define SANS_UMFPACK_TRANSPOSE        umfpack_di_transpose
#endif

} //namespace SLA
} //namespace numpack 

#endif //UMFPACKSOLVER_DEFINES_H
