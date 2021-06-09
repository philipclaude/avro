// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "numerics/matd.h"
#include "numerics/determinant.h"
#include "numerics/dual.h"
#include "numerics/symatd.h"
#include "numerics/tools.h"

#include "numerics/surreal/config.h"
#include "numerics/surreal/SurrealD.h"
#include "numerics/surreal/SurrealS.h"

#include "common/error.h"

#include <iostream>
#include <cmath>
#include <cstdio>
#include <vector>

#include <algorithm>
#include <numeric>

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

//extern "C" void dgetrs_( const char* TRANSA, int* N , double* A , int *LDA, int *IPIV , double *work , int* lwork , int *info );
//extern "C" void DGETRS( const char* TRANS, int* N , double* A , int *LDA, int *IPIV , double *work , int* lwork , int *info );

extern "C" void dgesvd_( char *jobu , char* jobvt , int *m , int *n , double *A , int* lda , double *s , double *u , int* ldu , double *vt , int* ldvt, double* work , int* lwork , int* info );

namespace avro
{

namespace numerics
{

template <typename type>
matd<type>::matd()
{
	m = 0;
	n = 0;
}

template <typename type>
matd<type>::matd(int _m,int _n)
{
	m = _m;
	n = _n;
	allocate();
}

template<typename type>
matd<type>::matd( const symatd<type>& T )
{
  m = T.n();
  n = m;

  allocate();

  for (int i=0;i<m;i++)
  for (int j=0;j<m;j++)
    operator()(i,j) = T(i,j);
}

template <typename type>
void
matd<type>::assertInRange(const int i,const int j)
{
	avro_assert_msg( i<m && j<n , "i = %d (m = %d) , j = %d (n = %d)",i,m,j,n);
}

template <typename type>
void
matd<type>::allocate()
{
	data.resize(m*n);
	zeros();
}

template <typename type>
void
matd<type>::zeros()
{
	int i,j;
	for (i=0;i<m;i++)
		for (j=0;j<n;j++)
			data[i*n+j] = 0.;
	return;
}

template <typename type>
void
matd<type>::eye()
{
	avro_assert_msg(m==n,"matrix not square: %dx%d",m,n);

	int i;
	zeros();
	for (i=0;i<m;i++)
		data[i*n+i] = type(1);

	return;
}

template <typename type>
void
matd<type>::solveLU(type *b)
{
	/*
		Solves (this)*x = b using Gauss elimination
		and back substitution with partial pivoting.
		The solution is overwritten into b. matd
		entries are also overwritten.
	*/
	int  i,j,ik,row,maxk;
	type factor,tmp=0,maxval;

	avro_assert_msg( m==n , "matrix not square: %dx%d",m,n);

	/* Gauss elimination */
	for (row=0;row<m;row++) {
		/* Find the row with the maximum value */
		maxval = data[row*n+row];
		maxk   = row;
		for (ik=row;ik<m;ik++) {
			if (fabs(data[ik*n+row]) > maxval) {
				maxk   = ik;
				maxval = fabs(data[ik*n+row]);
			}
		}

		/* Interchange the current row with the maxk'th row */
		for (ik=row;ik<m;ik++) {
			tmp = data[row*n+ik];
			data[row*n+ik]  = data[maxk*n+ik];
			data[maxk*n+ik] = tmp;
			tmp = b[row];
		}
		b[row]  = b[maxk];
		b[maxk] = tmp;

		factor = 1./data[row*n+row];
		for (j=0;j<m;j++) {
			data[row*n+j] = factor*data[row*n+j];
		}
		b[row] = factor*b[row];
		for (i=row+1;i<m;i++) {
			factor = data[i*n+row];
			for (j=0;j<m;j++) {
				data[i*n+j] = data[i*n+j] -factor*data[row*n+j];
			}
			b[i] = b[i] -factor*b[row];
		}
	}

	/* Back substitution */
	for (j=m-1;j>=1;j--) {
		for (i=j-1;i>=0;i--) {
			b[i] = b[i] -data[i*n+j]*b[j];
		}
	}

	/* Check for NaN */
	for (i=0;i<m;i++)
		avro_assert_msg( b[i]==b[i] , "nan :(" );
	return;
}

template <>
int
matd<double>::solve(double *x)
{
	avro_assert_msg( m==n , "should be square but is %dx%d",m,n );

  int nrhs = 1;
  int lda = m,ldb = m;
  int info;
	std::vector<int> ipiv(m);
	std::vector<double> tdata(m*m);

	for (int i=0;i<m;i++)
	for (int j=0;j<m;j++)
		tdata[i*m+j] = data[j*m+i];

#ifdef __APPLE__
  DGESV(&m,&nrhs,tdata.data(),&lda,ipiv.data(),x,&ldb,&info);
#else
  dgesv_(&m,&nrhs,tdata.data(),&lda,ipiv.data(),x,&ldb,&info);
#endif
	//if (info!=0) printf("info = %d\n",info);
	return info;
}

template<>
int
matd<double>::mrdivide( const matd<double>& Y , matd<double>& M )
{
	avro_assert_msg( m==n , "should be square but is %dx%d",m,n );

	int nrhs = Y.n;
	int lda = m,ldb = m;
	int info;
	std::vector<int> ipiv(m);
	std::vector<double> tdata(m*m);
	std::vector<double> mt(m*m);

	for (int i=0;i<m;i++)
	for (int j=0;j<m;j++)
	{
		tdata[i*m+j] = data[i*m+j];
		mt[i*m+j] = Y.data[i*m+j];
	}
#ifdef __APPLE__
	DGESV(&m,&nrhs,tdata.data(),&lda,ipiv.data(),mt.data(),&ldb,&info);
#else
	dgesv_(&m,&nrhs,tdata.data(),&lda,ipiv.data(),mt.data(),&ldb,&info);
#endif
	for (int i=0;i<m;i++)
	for (int j=0;j<m;j++)
		M.data[i*m+j] = mt[i*m+j];

	if (info!=0) printf("info = %d\n",info);
	return info;
}

template<>
int
matd<double>::kernel( matd<double>& K )
{
	int info;

	char jobu = 'A';
	char jobvt = 'A';

	int M = m;
	int N = n;

	int lda = m;
	int ldu = m;
	int ldvt = n;
	std::vector<double> S( m*n , 0. );
	std::vector<double> U( m*m , 0. );
	std::vector<double> VT( n*n , 0. );

	int lwork = std::max( 3*std::min(M,N)+std::max(M,N),5*std::min(M,N)-4 ) +10;
	std::vector<double> work( lwork );

	std::vector<double> tdata( m*n );
	for (int i=0;i<m;i++)
	for (int j=0;j<n;j++)
		tdata[j*m+i] = operator()(i,j);

  // perform the svd
	dgesvd_( &jobu , &jobvt , &M , &N , tdata.data() , &lda , S.data() , U.data() , &ldu , VT.data() , &ldvt , work.data() , &lwork , &info );

	// now we analyze the result to get the nullspace
	std::vector<double> s;
	if (m==0) s.push_back( S[0] );
	else if (m>0)
	{
		for (int i=0;i<m;i++)
			s.push_back( S[i] );
	}
	else s.push_back(0.);


	double maxS = *std::max_element( s.begin() , s.end() );
	double tol = std::max(m,n)*eps(maxS);
	int r = 0;
	for (int i=0;i<int(s.size());i++)
	{
		if (s[i]>tol) r++;
	}

	matd<double> v(N,N);
	for (int i=0;i<n;i++)
	for (int j=0;j<n;j++)
		v(i,j) = VT[j*n+i];

	// save the result
	K.m = n;
	K.n = n -r;
	K.allocate();
	for (int i=0;i<n;i++)
	for (int j=r;j<n;j++)
		K(i,j-r) = v(j,i);

	return info;
}

template<>
std::vector<double>
matd<double>::eig( matd<double>& q )
{
	avro_assert_msg( m==n , "should be square but is %dx%d",m,n );

	char jobz = 'V'; // N to get eigenvectors, "V" for no eigenvectors
	char uplo = 'L';
	int info;
	int lwork = m*(m+2);
	double norm;

	std::vector<double> work(m*(m+2));
	std::vector<double> lambda(m);
	std::vector<double> E(m*m);

	for (int i=0;i<m;i++)
	for (int j=0;j<m;j++)
		E[i*m+j] = data[j*m+i];

//#ifdef __APPLE__
#if 1
	dsyev_(&jobz,&uplo,&m,E.data(),&m,lambda.data(),work.data(),&lwork,&info);
#else
	DSYEV(&jobz,&uplo,&m,E.data(),&m,lambda.data(),work.data(),&lwork,&info);
#endif

	for (int j=0;j<m;j++)
	for (int i=0;i<m;i++)
		q(j,i) = E[i*m+j];

	if (info!=0) printf("info = %d\n",info);
	for (int i=0;i<m;i++)
	{
		norm = 0.;
		for (int j=0;j<m;j++)
			norm += q(j,i)*q(j,i);
		norm = sqrt(norm);
		for (int j=0;j<m;j++)
			q(j,i) /= norm;
	}
	return lambda;
}

template<>
void
matd<double>::multiply( const double* x , double *y ) const
{
	const char trans = 'T';
	const int inc = 1;
	const double a = 1.;
	const double b = 0.;

#ifdef __APPLE__
	DGEMV(&trans,&n,&m,&a,data.data(),&n,x,&inc,&b,y,&inc);
#else
	dgemv_(&trans,&n,&m,&a,data.data(),&n,x,&inc,&b,y,&inc);
#endif
}

template<>
void
matd<double>::multiply( const matd<double>& B , matd<double>& C ) const
{
	const char no_trans = 'N';
	const double a = 1.;
	const double b = 0.;

	int p = C.n;

	avro_assert( n==B.m );
	avro_assert( m==C.m );
	avro_assert( B.n==p );

	int m0 = m;
	int k0 = n;
	int n0 = B.n;

	try
	{
#ifdef __APPLE__
	DGEMM(&no_trans,&no_trans,&n0,&m0,&k0,&a,B.data.data(),&n0,data.data(),&k0,&b,C.data.data(),&n0);
#else
	dgemm_(&no_trans,&no_trans,&n0,&m0,&k0,&a,B.data.data(),&n0,data.data(),&k0,&b,C.data.data(),&n0);
#endif
	}
	catch(...)
	{
		printf("there was a dgemm error :(\n");
		avro_assert_not_reached;
	}
}

template<typename type>
void
matd<type>::multiply( const matd<type>& B , matd<type>& C ) const
{
  int p = C.n;
  avro_assert( n==B.m );
  avro_assert( m==C.m );
  avro_assert( B.n==p );

  for (int i=0;i<C.m;i++)
  for (int j=0;j<C.n;j++)
    C(i,j) = 0;

  for (int i=0;i<m;i++)
  for (int j=0;j<C.n;j++)
  for (int k=0;k<n;k++)
    C(i,j) += operator()(i,k)*B(k,j);
}

template<typename type>
matd<type>
matd<type>::operator*( const matd<type>& B ) const
{
	matd<type> C( m , B.n );
	multiply(B,C);
	return C;
}

template<typename type>
matd<type>
matd<type>::operator*( const symatd<type>& B0 ) const
{
	matd<type> B(B0);
	matd<type> C( m , B.n );
	multiply( B ,C);
	return C;
}

template <typename type>
void
matd<type>::solveLU(matd<type>& B)
{
	/*
		Solve Ax = B.
	*/
	int i,j,k;
	matd<type> A(m,n);
	std::vector<type> b(m);

	avro_assert( m==B.m && n==B.n );

	for (j=0;j<B.n;j++) {
		for (i=0;i<m;i++)
			for (k=0;k<n;k++)
				A(i,k) = operator()(i,k);

		for (i=0;i<m;i++) {
			b[i] = B(i,j);
		}
		A.solveLU(b.data());
		for (i=0;i<m;i++)
			B(i,j) = b[i];
	}
	return;
}

template<typename type>
void
matd<type>::inv( matd<type>& Minv ) const
{
	avro_assert( m==n );
	type idetM = 1./determinant(*this);

	if (n==1)
	{
		Minv(0,0) = idetM;
	}
	else if (n==2)
	{
		Minv(0,0) = operator()(1,1)*idetM;
		Minv(0,1) = -operator()(0,1)*idetM;
		Minv(1,0) = -operator()(1,0)*idetM;
		Minv(1,1) = operator()(0,0)*idetM;
	}
	else if (n==3)
	{
		type a1_1 = operator()(0,0); type a1_2 = operator()(0,1); type a1_3 = operator()(0,2);
		type a2_1 = operator()(1,0); type a2_2 = operator()(1,1); type a2_3 = operator()(1,2);
		type a3_1 = operator()(2,0); type a3_2 = operator()(2,1); type a3_3 = operator()(2,2);
		Minv(0,0) = (a2_2*a3_3 -a2_3*a3_2)*idetM;
		Minv(0,1) = (a1_3*a3_2 -a1_2*a3_3)*idetM;
		Minv(0,2) = (a1_2*a2_3 -a1_3*a2_2)*idetM;
		Minv(1,0) = (a2_3*a3_1 -a2_1*a3_3)*idetM;
		Minv(1,1) = (a1_1*a3_3 -a1_3*a3_1)*idetM;
		Minv(1,2) = (a1_3*a2_1 -a1_1*a2_3)*idetM;
		Minv(2,0) = (a2_1*a3_2 -a2_2*a3_1)*idetM;
		Minv(2,1) = (a1_2*a3_1 -a1_1*a3_2)*idetM;
		Minv(2,2) = (a1_1*a2_2 -a1_2*a2_1)*idetM;
	}
	else if (n==4)
	{
		type a1_1 = operator()(0,0); type a1_2 = operator()(0,1); type a1_3 = operator()(0,2); type a1_4 = operator()(0,3);
		type a2_1 = operator()(1,0); type a2_2 = operator()(1,1); type a2_3 = operator()(1,2); type a2_4 = operator()(1,3);
		type a3_1 = operator()(2,0); type a3_2 = operator()(2,1); type a3_3 = operator()(2,2); type a3_4 = operator()(2,3);
		type a4_1 = operator()(3,0); type a4_2 = operator()(3,1); type a4_3 = operator()(3,2); type a4_4 = operator()(3,3);

		Minv(0,0) = (a2_2*a3_3*a4_4-a2_2*a3_4*a4_3-a2_3*a3_2*a4_4+a2_3*a3_4*a4_2+a2_4*a3_2*a4_3-a2_4*a3_3*a4_2)*idetM;
		Minv(0,1) = (-a1_2*a3_3*a4_4+a1_2*a3_4*a4_3+a1_3*a3_2*a4_4-a1_3*a3_4*a4_2-a1_4*a3_2*a4_3+a1_4*a3_3*a4_2)*idetM;
		Minv(0,2) = (a1_2*a2_3*a4_4-a1_2*a2_4*a4_3-a1_3*a2_2*a4_4+a1_3*a2_4*a4_2+a1_4*a2_2*a4_3-a1_4*a2_3*a4_2)*idetM;
		Minv(0,3) = (-a1_2*a2_3*a3_4+a1_2*a2_4*a3_3+a1_3*a2_2*a3_4-a1_3*a2_4*a3_2-a1_4*a2_2*a3_3+a1_4*a2_3*a3_2)*idetM;
		Minv(1,0) = (-a2_1*a3_3*a4_4+a2_1*a3_4*a4_3+a2_3*a3_1*a4_4-a2_3*a3_4*a4_1-a2_4*a3_1*a4_3+a2_4*a3_3*a4_1)*idetM;
		Minv(1,1) = (a1_1*a3_3*a4_4-a1_1*a3_4*a4_3-a1_3*a3_1*a4_4+a1_3*a3_4*a4_1+a1_4*a3_1*a4_3-a1_4*a3_3*a4_1)*idetM;
		Minv(1,2) = (-a1_1*a2_3*a4_4+a1_1*a2_4*a4_3+a1_3*a2_1*a4_4-a1_3*a2_4*a4_1-a1_4*a2_1*a4_3+a1_4*a2_3*a4_1)*idetM;
		Minv(1,3) = (a1_1*a2_3*a3_4-a1_1*a2_4*a3_3-a1_3*a2_1*a3_4+a1_3*a2_4*a3_1+a1_4*a2_1*a3_3-a1_4*a2_3*a3_1)*idetM;
		Minv(2,0) = (a2_1*a3_2*a4_4-a2_1*a3_4*a4_2-a2_2*a3_1*a4_4+a2_2*a3_4*a4_1+a2_4*a3_1*a4_2-a2_4*a3_2*a4_1)*idetM;
		Minv(2,1) = (-a1_1*a3_2*a4_4+a1_1*a3_4*a4_2+a1_2*a3_1*a4_4-a1_2*a3_4*a4_1-a1_4*a3_1*a4_2+a1_4*a3_2*a4_1)*idetM;
		Minv(2,2) = (a1_1*a2_2*a4_4-a1_1*a2_4*a4_2-a1_2*a2_1*a4_4+a1_2*a2_4*a4_1+a1_4*a2_1*a4_2-a1_4*a2_2*a4_1)*idetM;
		Minv(2,3) = (-a1_1*a2_2*a3_4+a1_1*a2_4*a3_2+a1_2*a2_1*a3_4-a1_2*a2_4*a3_1-a1_4*a2_1*a3_2+a1_4*a2_2*a3_1)*idetM;
		Minv(3,0) = (-a2_1*a3_2*a4_3+a2_1*a3_3*a4_2+a2_2*a3_1*a4_3-a2_2*a3_3*a4_1-a2_3*a3_1*a4_2+a2_3*a3_2*a4_1)*idetM;
		Minv(3,1) = (a1_1*a3_2*a4_3-a1_1*a3_3*a4_2-a1_2*a3_1*a4_3+a1_2*a3_3*a4_1+a1_3*a3_1*a4_2-a1_3*a3_2*a4_1)*idetM;
		Minv(3,2) = (-a1_1*a2_2*a4_3+a1_1*a2_3*a4_2+a1_2*a2_1*a4_3-a1_2*a2_3*a4_1-a1_3*a2_1*a4_2+a1_3*a2_2*a4_1)*idetM;
		Minv(3,3) = (a1_1*a2_2*a3_3-a1_1*a2_3*a3_2-a1_2*a2_1*a3_3+a1_2*a2_3*a3_1+a1_3*a2_1*a3_2-a1_3*a2_2*a3_1)*idetM;
	}
	else
		avro_implement;
}

template<typename type>
void
matd<type>::column( index_t j , std::vector<type>& col ) const
{
  col.resize(m);
  for (int i=0;i<m;i++)
    col[i] = operator()(i,j);
}

template <typename type>
type&
matd<type>::operator() (const unsigned i,const unsigned j)
{
	assertInRange(i,j);
	return data[i*n+j];
}

template <typename type>
type
matd<type>::operator() (const unsigned i,const unsigned j) const
{
	return data[i*n+j];
}

template <typename type>
type&
matd<type>::val(const unsigned i,const unsigned j)
{
	assertInRange(i,j);
	return data[i*n+j];
}

template<typename type>
void
matd<type>::print()
{
	int i,j;
	for (i=0;i<m;i++) {
		printf("[ ");
		for (j=0;j<n;j++)
			std::cout << operator()(i,j) << " ";
		printf("]\n");
	}
	return;
}

template <typename type>
matd<type>::~matd()
{
	return;
}

template<typename type>
void
matd<type>::copy( const matd<type>& A , bool transpose )
{
	m = A.m;
	n = A.n;
	allocate();
	for (int i=0;i<m;i++)
	for (int j=0;j<n;j++)
	{
		if (transpose)
			operator()(j,i) = A(i,j);
		else
			operator()(i,j) = A(i,j);
	}
}

template<typename type>
matd<type>
matd<type>::transpose() const
{
	matd<type> B(n,m);
	for (int i=0;i<m;i++)
	for (int j=0;j<n;j++)
		B(j,i) = operator()(i,j);
	return B;
}

template<typename type>
matd<type>
matd<type>::inverse() const
{
	avro_assert( m==n );
	matd<type> Ainv(m,n);
	inv(Ainv);
	return Ainv;
}

template <>
void
matd<double>::display() const
{
	int i,j;
	for (i=0;i<m;i++) {
		for (j=0;j<n;j++) {
			printf("%3.6E ",(double)data[i*n+j]);
		}
		printf("\n");
	}
	return;
}

template<>
void
matd<double>::tomatlab( const std::string& var ) const
{
  printf("%s = [",var.c_str());
  for (int i=0;i<m;i++)
  {
    for (int j=0;j<n;j++)
    {
      printf("%1.16e",operator()(i,j));
      if (j<n-1) printf(",");
    }
    if (i<m-1)
      printf(";");
  }
  printf("];\n");
}

template<typename type>
type
matd<type>::det()
{
  // get the eigenvalues and multiple them
  avro_implement;
}

template<typename type>
void
matd<type>::rankone( std::vector<type>& x )
{
  avro_assert( m==n );
  avro_assert( x.size()==index_t(m) );

  for (int i=0;i<m;i++)
  for (int j=0;j<m;j++)
    operator()(i,j) = x[i]*x[j];
}

template<typename type>
void
matd<type>::scale( const type& alpha )
{
  for (index_t k=0;k<data.size();k++)
    data[k] *= alpha;
}

template<typename type>
void
matd<type>::add( const matd<type>& x )
{
  avro_assert( x.m==m );
  avro_assert( x.n==n );

  for (int i=0;i<m;i++)
  for (int j=0;j<n;j++)
    operator()(i,j) += x(i,j);
}

template<>
void
matd<double>::range( matd<double>& U0 ) const
{
	int info;

	char jobu = 'A';
	char jobvt = 'A';

	int M = m;
	int N = n;

	int lda = m;
	int ldu = m;
	int ldvt = n;
	std::vector<double> S( m*n , 0. );
	std::vector<double> U( m*m , 0. );
	std::vector<double> VT( n*n , 0. );

	int lwork = std::max( 3*std::min(M,N)+std::max(M,N),5*std::min(M,N)-4 ) +10;
	std::vector<double> work( lwork );

	std::vector<double> tdata( m*n );
	for (int i=0;i<m;i++)
	for (int j=0;j<n;j++)
		tdata[j*m+i] = operator()(i,j);

	// perform the svd
	dgesvd_( &jobu , &jobvt , &M , &N , tdata.data() , &lda , S.data() , U.data() , &ldu , VT.data() , &ldvt , work.data() , &lwork , &info );

	// analyze the singular values to determine the rank
	std::vector<double> s;
	if (m==0) s.push_back( S[0] );
	else if (m>0)
	{
		for (int i=0;i<m;i++)
			s.push_back( S[i] );
	}
	else s.push_back(0.);

	// compute the tolerance
	double maxS = *std::max_element( s.begin() , s.end() );
	double tol = std::max(m,n)*eps(maxS);
	int r = 0;
	for (int i=0;i<int(s.size());i++)
	{
		if (s[i]>tol) r++;
	}

	U0.m = M;
	U0.n = r;
	U0.allocate();

	// save the left singular vectors which have non-zero singular values
	index_t k = 0;
	for (index_t j=0;j<index_t(M);j++)
	{
		if (s[j]<=tol) continue;
		//printf("left singular vector for svalue %g\n",s[j]);
		for (index_t i=0;i<index_t(M);i++)
		{
			U0(i,k) = U[j*m+i];
		}
		k++;
	}
}

template class matd<double>;
template class matd<dual>;

#if USE_SURREAL
template class matd<SurrealS<1>>;
template class matd<SurrealS<3>>;
template class matd<SurrealS<6>>;
template class matd<SurrealS<10>>;
#endif

} // numerics

} // avro
