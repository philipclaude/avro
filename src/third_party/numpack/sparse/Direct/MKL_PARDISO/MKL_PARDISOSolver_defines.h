// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MKL_PARDISOSOLVER_DEFINES_H
#define MKL_PARDISOSOLVER_DEFINES_H

// MKL_PARDISO interface defines

#ifdef INTEL_MKL

#include <mkl_types.h>

#define USE_SANS_MKL_PARDISO_LONG

#ifdef USE_SANS_MKL_PARDISO_LONG            // MKL_PARDISO 64-bit interface
  #define SANS_MKL_PARDISO_INT              MKL_INT64
  #define SANS_MKL_PARDISO                  PARDISO_64

#else                                       // MKL_PARDISO 32-bit interface
  #define SANS_MKL_PARDISO_INT              MKL_INT
  #define SANS_MKL_PARDISO                  PARDISO
#endif

#else //INTEL_MKL

#define SANS_MKL_PARDISO_INT long

#endif //INTEL_MKL

#endif //MKL_PARDISOSOLVER_DEFINES_H
