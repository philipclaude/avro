//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_NUMERICS_LAPACK_H_
#define avro_LIB_NUMERICS_LAPACK_H_

namespace avro
{

namespace numerics
{

// solving linear system
extern "C" void DGESV(int *N, int *NRHS, double *A, int *LDA, int *IPIV,double *B, int *LDB, int *INFO);
extern "C" void dgesv_(int *N, int *NRHS, double *A, int *LDA, int *IPIV,double *B, int *LDB, int *INFO);

// eigenvalues
extern "C" void DSYEV(char *jobz,char *uplo,int *m, double *E, int *lda,double *w,double *work,int* lwork,int *info);
extern "C" void dsyev_(char *jobz,char *uplo,int *m, double *E, int *lda,double *w,double *work,int* lwork,int *info);

// multiply matrix by vector
extern "C" void dgemv_(const char *TRANS, const int* M,const int *N,const double* ALPHA,const double* A, const int *LDA,const double *X,const int* INCX,const double* BETA,double *Y,const int *INCY);
extern "C" void DGEMV(const char *TRANS, const int* M,const int *N,const double* ALPHA,const double* A, const int *LDA,const double *X,const int* INCX,const double* BETA,double *Y,const int *INCY);

// multiply matrix by matrix
extern "C" void dgemm_(const char *TRANSA,const char* TRANSB, const int *M, const int *N,const int *K , const double* alpha, const double *A,const int *lda,const double* B, const int *ldb, const double *beta, double *C, const int *ldc);
extern "C" void DGEMM(const char *TRANSA,const char* TRANSB, const int *M, const int *N,const int *K , const double* alpha, const double *A,const int *lda,const double* B, const int *ldb, const double *beta, double *C, const int *ldc);

// LU factorization
extern "C" void dgetrf_( int* M , int* N , double *A , int* LDA , int* IPIV , int* INFO );
extern "C" void DGETRF( int* M , int* N , double *A , int* LDA , int* IPIV , int* INFO );

// SVD
extern "C" void dgesvd_( char *jobu , char* jobvt , int *m , int *n , double *A , int* lda , double *s , double *u , int* ldu , double *vt , int* ldvt, double* work , int* lwork , int* info );
extern "C" void DGESVD( char *jobu , char* jobvt , int *m , int *n , double *A , int* lda , double *s , double *u , int* ldu , double *vt , int* ldvt, double* work , int* lwork , int* info );

// SVD
extern "C" void dgeev_( char* jobvl, char* jobvr, int* n, double* a,
                   int* lda, double* wr, double* wi, double* vl,
                   int* ldvl, double* vr, int* ldvr, double* work,
                   int* lwork, int *info );
extern "C" void DGEEV( char* jobvl, char* jobvr, int* n, double* a,
                   int* lda, double* wr, double* wi, double* vl,
                   int* ldvl, double* vr, int* ldvr, double* work,
                   int* lwork, int *info );


} // numerics

} // avro

#ifdef AVRO_NO_LAPACK

extern "C"
{

void
dgesv_(int *N, int *NRHS, double *A, int *LDA,
       int *IPIV,double *B, int *LDB, int *INFO) {
  avro_assert_not_reached;
}

void dgesvd_( char *jobu , char* jobvt , int *m , int *n , double* a ,
              int* lda , double *s , double *u , int* ldu , double *vt ,
              int* ldvt, double* work , int* lwork , int* info ) {
  avro_assert_not_reached;
}

void
dgeev_( char* jobvl, char* jobvr, int* n, double* a,
        int* lda, double* wr, double* wi, double* vl,
        int* ldvl, double* vr, int* ldvr, double* work,
        int* lwork, int *info ) {
  avro_assert_not_reached;
}

}
#endif

#endif
