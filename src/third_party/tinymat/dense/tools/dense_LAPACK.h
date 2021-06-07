// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef dense_LAPACK_H
#define dense_LAPACK_H

// Note that none of the LAPACKE functions should be used other than with testing.
// LAPACKE causes problems linking with the OS X Accelerate library.

#define DLA_LAPACK
#ifdef DLA_LAPACK

#ifdef LAPACK_ACCELERATE

#include <Accelerate/Accelerate.h>

#define lapack_int __CLPK_integer
#define lapack_logical __CLPK_logical

#define LAPACK_dgeev dgeev_
#define LAPACK_sgeev sgeev_

#define LAPACK_dgesvd dgesvd_
#define LAPACK_sgesvd sgesvd_

#elif defined(DLA_LAPACKE)

// Need these definitions because lapacke is not c++11 compatible yet. These #defs  gets around that for now..
#define lapack_complex_float float _Complex
#define lapack_complex_double double _Complex
#define LAPACK_COMPLEX_CPP
#include <lapacke.h>

#else

// Borrowed from lapacke logic lapacke_mangling.h
#ifndef LAPACK_GLOBAL
#if defined(LAPACK_GLOBAL_PATTERN_LC)
#define LAPACK_GLOBAL(lcname,UCNAME)  lcname##_
#elif defined(LAPACK_GLOBAL_PATTERN_UC)
#define LAPACK_GLOBAL(lcname,UCNAME)  UCNAME
#elif defined(LAPACK_GLOBAL_PATTERN_MC)
#define LAPACK_GLOBAL(lcname,UCNAME)  lcname
#else
#define LAPACK_GLOBAL(lcname,UCNAME)  lcname##_
#endif
#endif

// borrowed from lapacke.h
#ifndef lapack_int
#define lapack_int     int
#endif

extern "C"
{

#define LAPACK_dgeev LAPACK_GLOBAL(dgeev,DGEEV)
#define LAPACK_sgeev LAPACK_GLOBAL(sgeev,SGEEV)

void LAPACK_dgeev( char* jobvl, char* jobvr, lapack_int* n, double* a,
                   lapack_int* lda, double* wr, double* wi, double* vl,
                   lapack_int* ldvl, double* vr, lapack_int* ldvr, double* work,
                   lapack_int* lwork, lapack_int *info );
void LAPACK_sgeev( char* jobvl, char* jobvr, lapack_int* n, float* a,
                   lapack_int* lda, float* wr, float* wi, float* vl,
                   lapack_int* ldvl, float* vr, lapack_int* ldvr, float* work,
                   lapack_int* lwork, lapack_int *info );

#define LAPACK_dgesvd LAPACK_GLOBAL(dgesvd,DGESVD)
#define LAPACK_sgesvd LAPACK_GLOBAL(sgesvd,SGESVD)

void LAPACK_dgesvd( char* jobu, char* jobvt, lapack_int* m, lapack_int* n,
                    double* a, lapack_int* lda, double* s, double* u,
                    lapack_int* ldu, double* vt, lapack_int* ldvt, double* work,
                    lapack_int* lwork, lapack_int *info );
void LAPACK_sgesvd( char* jobu, char* jobvt, lapack_int* m, lapack_int* n,
                    float* a, lapack_int* lda, float* s, float* u,
                    lapack_int* ldu, float* vt, lapack_int* ldvt, float* work,
                    lapack_int* lwork, lapack_int *info );

} // extern "C"

#endif //LAPACK_ACCELERATE or DLA_LAPACKE

#endif //DLA_LAPACK

#endif //dense_LAPACK_H
