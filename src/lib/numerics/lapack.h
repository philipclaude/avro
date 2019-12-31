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

} // numerics

} // avro

#endif
