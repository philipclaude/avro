// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef DENSELINALG_BLAS_H
#define DENSELINALG_BLAS_H

#ifdef DLA_BLAS_ATLAS

#define DLA_BLAS
extern "C"
{
#include <cblas.h>
}

#elif defined(DLA_BLAS_GOTO)

#define DLA_BLAS
extern "C"
{
#include <common.h>
#include <cblas.h>
}

#elif defined(DLA_BLAS_MKL)

#define DLA_BLAS
#include <mkl_cblas.h>

#elif defined(DLA_BLAS_ACCELERATE)

#define DLA_BLAS
#include <Accelerate/Accelerate.h>

#elif defined(DLA_BLAS_GENERIC)

#define DLA_BLAS
extern "C"
{
#include <cblas.h>
}

#endif

#endif //DENSELINALG_BLAS_H
