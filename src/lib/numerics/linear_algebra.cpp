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
kernel( const MatrixD<real_t>& A , MatrixD<real_t>& K )
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

	MatrixD<real_t> v(N,N);
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
range( const MatrixD<real_t>& A , MatrixD<real_t>& U0 )
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
  #ifndef AVRO_NO_LAPACK
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

} // numerics

} // avro
