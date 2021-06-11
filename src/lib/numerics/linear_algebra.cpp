//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "common/tools.h"

#include "numerics/lapack.h"
#include "numerics/linear_algebra.h"
#include "numerics/mat.h"
#include "numerics/vec.h"

#include <algorithm>
#include <vector>

#ifdef AVRO_NO_LAPACK
extern "C"
{
void dgesv_(int *N, int *NRHS, double *A, int *LDA, int *IPIV,double *B, int *LDB, int *INFO)
{
avro_assert_not_reached;
}

extern "C" void dgesvd_( char *jobu , char* jobvt , int *m , int *n , double *A , int* lda , double *s , double *u , int* ldu , double *vt , int* ldvt, double* work , int* lwork , int* info )
{
avro_assert_not_reached;
}

extern "C"
void dgeev_( char* jobvl, char* jobvr, int* n, double* a,
                   int* lda, double* wr, double* wi, double* vl,
                   int* ldvl, double* vr, int* ldvr, double* work,
                   int* lwork, int *info )
{
  avro_assert_not_reached;
}

}
#endif

namespace avro
{

namespace numerics
{

real_t
eps( const real_t& x )
{
  // this is like matlab's eps function
  return ::nextafter( x , x +1.0f ) -x;
}

template<>
int
kernel( const matd<real_t>& A , matd<real_t>& K )
{
  int info;

	char jobu = 'A';
	char jobvt = 'A';

  int m = A.m();
  int n = A.n();

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
		tdata[j*m+i] = A(i,j);

  // perform the svd
        #ifdef AVRO_NO_LAPACK
        avro_assert_not_reached;
        #else
	dgesvd_( &jobu , &jobvt , &M , &N , tdata.data() , &lda , S.data() , U.data() , &ldu , VT.data() , &ldvt , work.data() , &lwork , &info );
	#endif

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

	matd<real_t> v(N,N);
	for (int i=0;i<n;i++)
	for (int j=0;j<n;j++)
		v(i,j) = VT[j*n+i];

	// save the result
	//K.m = n;
	//K.n = n -r;
	//K.allocate();
  K.resize(n,n-r);
	for (int i=0;i<n;i++)
	for (int j=r;j<n;j++)
		K(i,j-r) = v(j,i);

	return info;
}

template<>
int
range( const matd<real_t>& A , matd<real_t>& U0 )
{
  int info;

  char jobu = 'A';
  char jobvt = 'A';

  int M = A.m();
  int N = A.n();

  int m = M;
  int n = N;

  int lda = M;
  int ldu = M;
  int ldvt = N;
  std::vector<double> S( m*n , 0. );
  std::vector<double> U( m*m , 0. );
  std::vector<double> VT( n*n , 0. );

  int lwork = std::max( 3*std::min(M,N)+std::max(M,N),5*std::min(M,N)-4 ) +10;
  std::vector<double> work( lwork );

  std::vector<double> tdata( m*n );
  for (int i=0;i<m;i++)
  for (int j=0;j<n;j++)
    tdata[j*m+i] = A(i,j);

  // perform the svd
  #ifdef AVRO_NO_LAPACK
  UNUSED(jobu);
  UNUSED(jobvt);
  UNUSED(lda);
  UNUSED(ldu);
  UNUSED(ldvt);
  avro_assert_not_reached;
  #else
  dgesvd_( &jobu , &jobvt , &M , &N , tdata.data() , &lda , S.data() , U.data() , &ldu , VT.data() , &ldvt , work.data() , &lwork , &info );
  #endif

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

  //U0.m = M;
  //U0.n = r;
  //U0.allocate();
  U0.resize(M,r);

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
  return info;
}

template<typename T>
matd<T>
transpose( const matd<T>& A ) {
  matd<T> At( A.n() , A.m() );
  for (index_t i = 0; i < A.m(); i++)
  for (index_t j = 0; j < A.n(); j++)
    At(j,i) = A(i,j);
  return At;
}

template<typename T>
matd<T>
diag( const vecd<T>& d ) {
  matd<T> A( d.m() , d.m() ); // initializes to zero
  for (index_t i = 0; i < d.m(); i++)
    A(i,i) = d(i);
  return A;
}

template <typename T>
void
decomposeLUP(const matd<T>& A , matd<T>& LU , std::vector<index_t>& P)
{
  real_t tol = 1e-12;

  index_t m = A.m();
  index_t n = A.n();
	avro_assert_msg( m==n , "matrix not square: %lu x %lu",m,n);

  // define unit permutation matrix
  P.resize( n + 1 );
  for (index_t j = 0; j <= n; j++)
    P[j] = j;

  // copy A into B
  for (index_t i = 0; i < m; i++)
  for (index_t j = 0; j < n; j++)
    LU(i,j) = A(i,j);

  // initialize the rows that will be used to swap when pivoting
  vecd<T> row1(n);
  vecd<T> row2(n);
	for (index_t i = 0; i < m; i++) {

		// find the row with the maximum value
		T maxa       = 0;
		index_t imax = i;

    for (index_t k = i; k < n; k++) {
      T absa = fabs(LU(k,i));
      if (absa > maxa) {
        maxa = absa;
        imax = k;
      }
    }

    avro_assert_msg( maxa > tol , "matrix is degenerate" );

    if (imax != i) {
      // pivot P
      index_t j = P[i];
      P[i]      = P[imax];
      P[imax]   = j;

      // pivot rows of LU
      LU.get_row( i , row1 );
      LU.get_row( imax , row2 );
      LU.set_row( i , row2 );
      LU.set_row( imax , row1 );

      // increase the pivot counter for determinant
      P[n]++;
    }

    for (index_t j = i+1; j < n; j++) {
      LU(j,i) /= LU(i,i);

      for (index_t k = i+1; k < n; k++)
        LU(j,k) -= LU(j,i) * LU(i,k);
    }
  }
}

template<typename T>
void
solveLUP( const matd<T>& LU , const std::vector<index_t>& P , const vecd<T>& b , vecd<T>& x )
{
  index_t n = LU.n();
  avro_assert( x.m() == b.m() );
  avro_assert( n == LU.m() );

  // pivot the vector according to P
  for (index_t i = 0; i < n; i++) {
    x(i) = b( P[i] );
    for (index_t k = 0; k < i; k++)
      x(i) -= LU(i,k)*x(k);
  }

  // backwards substitution
  for (int i = n-1; i >= 0; i-- ) {
    for (index_t k = i+1; k < n; k++)
      x(i) -= LU(i,k)*x(k);
    x(i) /= LU(i,i);
  }
}

template<typename T>
void
solveLUP( const matd<T>& A , const vecd<T>& b , vecd<T>& x )
{
  matd<T> LU( A.m() , A.n() );
  std::vector<index_t> P(A.n()+1);
  decomposeLUP(A,LU,P);
  solveLUP(LU,P,b,x);
}

template<typename T>
void
inverseLUP( const matd<T>& LU , const std::vector<index_t>& P , matd<real_t>& Ainv )
{
  index_t n = LU.n();
  for (index_t j = 0; j < n; j++)
  {
    for (index_t i = 0; i < n; i++) {
      Ainv(i,j) = (P[i] == j) ? 1.0 : 0.0;

      for (index_t k = 0; k < i; k++)
        Ainv(i,j) -= LU(i,k) * Ainv(k,j);
    }

    for (index_t i = n-1; i >= 0; i--) {
      for (index_t k = i+1; k < n; k++)
        Ainv(i,j) -= LU(i,k)*Ainv(k,j);
      Ainv(i,j) /= LU(i,i);
    }
  }
}

template<typename T>
void
inverseLUP( const matd<T>& A , matd<T>& Ainv )
{
  matd<T> LU( A.m() , A.n() );
  std::vector<index_t> P(A.n()+1);
  decomposeLUP(A,LU,P);
  inverseLUP(LU,P,Ainv);
}

template void solveLUP( const matd<real_t>& , const vecd<real_t>& , vecd<real_t>& );
template void inverseLUP( const matd<real_t>& , matd<real_t>& );

#define INSTANTIATE_TRANSPOSE(T) template matd<T> numerics::transpose( const matd<T>& );
INSTANTIATE_TRANSPOSE( real_t )
INSTANTIATE_TRANSPOSE( SurrealS<1> )
INSTANTIATE_TRANSPOSE( SurrealS<3> )
INSTANTIATE_TRANSPOSE( SurrealS<6> )
INSTANTIATE_TRANSPOSE( SurrealS<10> )
#undef INSTANTIATE_TRANSPOSE

#define INSTANTIATE_DIAG(T) template matd<T> diag( const vecd<T>& );
INSTANTIATE_DIAG( real_t )
INSTANTIATE_DIAG(  SurrealS<1> )
INSTANTIATE_DIAG(  SurrealS<3> )
INSTANTIATE_DIAG(  SurrealS<6> )
INSTANTIATE_DIAG(  SurrealS<10> )
#undef INSTANTIATE_DIAG

} // numerics

} // avro
