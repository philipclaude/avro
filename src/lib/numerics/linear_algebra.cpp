//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "common/error.h"
#include "common/tools.h"
#include "types.h"

#include "numerics/mat.h"
#include "numerics/sym.h"
#include "numerics/dual.h"

#include "numerics/functions.h"
#include "numerics/lapack.h"
#include "numerics/linear_algebra.h"
#include "numerics/mat.h"
#include "numerics/vec.h"

#include "numerics/surreal/config.h"
#include "numerics/surreal/SurrealD.h"
#include "numerics/surreal/SurrealS.h"

#include <algorithm>
#include <vector>

namespace avro
{

namespace numerics
{

template<>
int
kernel( const matd<real_t>& A , matd<real_t>& K ) {
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
	for (int i = 0; i < m; i++)
	for (int j = 0; j < n; j++)
		tdata[j*m+i] = A(i,j);

  // perform the svd
  #ifdef AVRO_NO_LAPACK
  avro_assert_not_reached;
  #else
	dgesvd_( &jobu , &jobvt , &M , &N , tdata.data() , &lda , S.data() , U.data() , &ldu , VT.data() , &ldvt , work.data() , &lwork , &info );
	#endif

	// now we analyze the result to get the nullspace
	std::vector<double> s;
	if (m == 0) s.push_back( S[0] );
	else if (m > 0) {
		for (int i = 0; i < m; i++)
			s.push_back( S[i] );
	}
	else s.push_back(0.);


	double maxS = *std::max_element( s.begin() , s.end() );
	double tol = std::max(m,n)*eps(maxS);
	int r = 0;
	for (int i = 0; i < int(s.size()); i++) {
		if (s[i]>tol) r++;
	}

	matd<real_t> v(N,N);
	for (int i = 0; i < n;i++)
	for (int j = 0; j < n;j++)
		v(i,j) = VT[j*n+i];

	// save the result
  K.resize(n,n-r);
	for (int i = 0; i < n; i++)
	for (int j = r; j < n; j++)
		K(i,j-r) = v(j,i);

	return info;
}

template<>
int
range( const matd<real_t>& A , matd<real_t>& U0 ) {
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

  // TODO matd's are now column major, so we don't need temporary anymore
  std::vector<double> tdata( m*n );
  for (int i = 0; i < m; i++)
  for (int j = 0; j < n; j++)
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
  if (m == 0 ) s.push_back( S[0] );
  else if (m > 0)
  {
    for (int i =0 ; i < m; i++)
      s.push_back( S[i] );
  }
  else s.push_back(0.);

  // compute the tolerance
  double maxS = *std::max_element( s.begin() , s.end() );
  double tol = std::max(m,n)*eps(maxS);
  int r = 0;
  for (int i = 0; i < int(s.size()); i++) {
    if (s[i]>tol) r++;
  }
  U0.resize(M,r);

  // save the left singular vectors which have non-zero singular values
  index_t k = 0;
  for (index_t j = 0; j < index_t(M); j++)
  {
    if (s[j] <= tol) continue;
    //printf("left singular vector for svalue %g\n",s[j]);
    for (index_t i = 0; i < index_t(M); i++) {
      U0(i,k) = U[j*m+i];
    }
    k++;
  }
  return info;
}

template<typename T>
symd<T>
interp( const std::vector<real_t>& alpha , const std::vector<symd<T>>& tensors ) {
  avro_assert( alpha.size() == tensors.size() );
  avro_assert( tensors.size() > 0 );
  const index_t n = tensors[0].n();
  symd<T> m(n);
  m.zero();
  for (index_t k = 0; k < tensors.size(); k++)
    m = m +numerics::logm(tensors[k])*alpha[k];
  return numerics::expm(m);
}

template<typename T>
void
eign( const matd<T>& A , vecd<T>& L , matd<T>& Q ) {
  avro_implement;
}

template<typename T>
void
eign( const matd<T>& A , vecd<T>& L ) {
  avro_implement;
}

template<>
void
eign( const matd<real_t>& A , vecd<real_t>& L ) {

  avro_assert( A.m() == A.n()  );
  avro_assert( A.m() == L.m() );

  #if 0 //ndef AVRO_NO_LAPACK
  avro_implement;
  #else
  char jobvl = 'N'; // left eigenvector
  char jobvr = 'N'; // right eigenvector
  real_t *vl = NULL;     // left eigenvector place holder
  int ldvl = 1;
  real_t *vr = NULL;     // right eigenvector place holder
  int ldvr = 1;
  int m = A.m();
  int stride = m;
  int lwork = 16*A.m();
  std::vector<real_t> work( lwork );
  std::vector<real_t> wi( m );
  int INFO;

  dgeev_(&jobvl,&jobvr,&m,const_cast<real_t*>(&A(0,0)),&stride,&L[0],&wi[0],vl,&ldvl,vr,&ldvr,&work[0],&lwork,&INFO);
  avro_assert_msg( INFO == 0, "INFO == %d", INFO );
  #endif
}

template<>
std::pair< vecd<real_t> , matd<real_t> >
eign( const symd<real_t>& A ) {

  avro_assert( A.m() == A.n()  );
  int n = A.n();

  vecd<real_t> L(n);
  matd<real_t> Q(n,n);

  #if 0 //ndef AVRO_NO_LAPACK
  avro_implement;
  #else
  char jobvl = 'N';  // do not compute left eigenvectors
  char jobvr = 'V';  // compute right eigenvectors
  real_t *vl = NULL; // left eigenvector place holder
  int ldvl = 1;
  std::vector<real_t> vr(n*n);    // right eigenvector place holder
  int ldvr = n;
  int m = A.m();
  int stride = m;
  int lwork = 16*A.m();
  std::vector<real_t> work( lwork );
  std::vector<real_t> wi( m );
  int INFO;

  matd<real_t> Atmp(A);
  dgeev_(&jobvl,&jobvr,&m,&Atmp(0,0),&stride,&L[0],&wi[0],vl,&ldvl,vr.data(),&ldvr,&work[0],&lwork,&INFO);
  avro_assert_msg( INFO == 0, "INFO == %d", INFO );

  index_t k = 0;
  for (index_t j = 0; j < n; j++) // through columns
  for (index_t i = 0; i < n; i++) // through rows
    Q(i,j) = vr[k++];
  #endif
  return {L,Q};
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
solveLUP( const matd<T>& A , const vecd<T>& b , vecd<T>& x ) {
  matd<T> LU( A.m() , A.n() );
  std::vector<index_t> P(A.n()+1);
  decomposeLUP(A,LU,P);
  solveLUP(LU,P,b,x);
}

template<typename T>
void
inverseLUP( const matd<T>& LU , const std::vector<index_t>& P , matd<real_t>& Ainv ) {
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
matd<T>
inverse( const matd<T>& M )
{
	avro_assert( M.m() == M.n() );
	const T idetM = 1./det(M);

  matd<T> Minv(M.m(),M.n());

	if (M.n() == 1) {
		Minv(0,0) = idetM;
	}
	else if (M.n() == 2) {
		Minv(0,0) =  M(1,1)*idetM;
		Minv(0,1) = -M(0,1)*idetM;
		Minv(1,0) = -M(1,0)*idetM;
		Minv(1,1) =  M(0,0)*idetM;
	}
	else if (M.n() == 3) {
		const T a1_1 = M(0,0); const T a1_2 = M(0,1); const T a1_3 = M(0,2);
		const T a2_1 = M(1,0); const T a2_2 = M(1,1); const T a2_3 = M(1,2);
		const T a3_1 = M(2,0); const T a3_2 = M(2,1); const T a3_3 = M(2,2);
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
	else if (M.n() == 4) {
		const T a1_1 = M(0,0); const T a1_2 = M(0,1); const T a1_3 = M(0,2); const T a1_4 = M(0,3);
		const T a2_1 = M(1,0); const T a2_2 = M(1,1); const T a2_3 = M(1,2); const T a2_4 = M(1,3);
		const T a3_1 = M(2,0); const T a3_2 = M(2,1); const T a3_3 = M(2,2); const T a3_4 = M(2,3);
		const T a4_1 = M(3,0); const T a4_2 = M(3,1); const T a4_3 = M(3,2); const T a4_4 = M(3,3);

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
  return Minv;
}

template<typename T>
symd<T>
inverse( const symd<T>& M )
{
	const T idetM = 1./det(M);

  symd<T> Minv(M.m(),M.n());

	if (M.n() == 1)
	{
		Minv(0,0) = idetM;
	}
	else if (M.n() == 2)
	{
		Minv(0,0) =  M(1,1)*idetM;
		Minv(0,1) = -M(0,1)*idetM;
		Minv(1,0) = -M(1,0)*idetM;
		Minv(1,1) =  M(0,0)*idetM;
	}
	else if (M.n() == 3)
	{
		const T a1_1 = M(0,0); const T a1_2 = M(0,1); const T a1_3 = M(0,2);
		const T a2_1 = M(1,0); const T a2_2 = M(1,1); const T a2_3 = M(1,2);
		const T a3_1 = M(2,0); const T a3_2 = M(2,1); const T a3_3 = M(2,2);
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
	else if (M.n()==4)
	{
		const T a1_1 = M(0,0); const T a1_2 = M(0,1); const T a1_3 = M(0,2); const T a1_4 = M(0,3);
		const T a2_1 = M(1,0); const T a2_2 = M(1,1); const T a2_3 = M(1,2); const T a2_4 = M(1,3);
		const T a3_1 = M(2,0); const T a3_2 = M(2,1); const T a3_3 = M(2,2); const T a3_4 = M(2,3);
		const T a4_1 = M(3,0); const T a4_2 = M(3,1); const T a4_3 = M(3,2); const T a4_4 = M(3,3);

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
  return Minv;
}

template<typename T>
std::pair< vecd<T> , matd<T> >
eig( const symd<T>& m ) {
  return m.eig();
}

template<typename T>
void
eig( const symd<T>& m , vecd<T>& L , matd<T>& Q ) {
  std::pair< vecd<T> , matd<T> > decomp = m.eig();
  L.set( decomp.first );
  Q.set( decomp.second );
}

template<typename T>
void
inverseLUP( const matd<T>& A , matd<T>& Ainv ) {
  matd<T> LU( A.m() , A.n() );
  std::vector<index_t> P(A.n()+1);
  decomposeLUP(A,LU,P);
  inverseLUP(LU,P,Ainv);
}


template<typename type>
type
det(const matd<type>& X)
{
  avro_assert(X.m()==X.n());
  if (X.m()==1) return X(0,0);
  if (X.m()==2) return X(1,1)*X(0,0)-X(0,1)*X(1,0);
  if (X.m()==3) return X(0,0)*(X(2,2)*X(1,1)-X(2,1)*X(1,2))
                    -X(1,0)*(X(2,2)*X(0,1)-X(2,1)*X(0,2))
                    +X(2,0)*(X(1,2)*X(0,1)-X(1,1)*X(0,2));
  const type X1_1=X(0,0),X1_2=X(0,1),X1_3=X(0,2),X1_4=X(0,3);
  const type X2_1=X(1,0),X2_2=X(1,1),X2_3=X(1,2),X2_4=X(1,3);
  const type X3_1=X(2,0),X3_2=X(2,1),X3_3=X(2,2),X3_4=X(2,3);
  const type X4_1=X(3,0),X4_2=X(3,1),X4_3=X(3,2),X4_4=X(3,3);
  if (X.m()==4)
  {   return
        X1_1*X2_2*X3_3*X4_4 - X1_1*X2_2*X3_4*X4_3 - X1_1*X2_3*X3_2*X4_4
      + X1_1*X2_3*X3_4*X4_2 + X1_1*X2_4*X3_2*X4_3 - X1_1*X2_4*X3_3*X4_2
      - X1_2*X2_1*X3_3*X4_4 + X1_2*X2_1*X3_4*X4_3 + X1_2*X2_3*X3_1*X4_4
      - X1_2*X2_3*X3_4*X4_1 - X1_2*X2_4*X3_1*X4_3 + X1_2*X2_4*X3_3*X4_1
      + X1_3*X2_1*X3_2*X4_4 - X1_3*X2_1*X3_4*X4_2 - X1_3*X2_2*X3_1*X4_4
      + X1_3*X2_2*X3_4*X4_1 + X1_3*X2_4*X3_1*X4_2 - X1_3*X2_4*X3_2*X4_1
      - X1_4*X2_1*X3_2*X4_3 + X1_4*X2_1*X3_3*X4_2 + X1_4*X2_2*X3_1*X4_3
      - X1_4*X2_2*X3_3*X4_1 - X1_4*X2_3*X3_1*X4_2 + X1_4*X2_3*X3_2*X4_1;
  }
  const type X1_5=X(0,4),X2_5=X(1,4),X3_5=X(2,4),X4_5=X(3,4);
  const type X5_1=X(4,0),X5_2=X(4,1),X5_3=X(4,2),X5_4=X(4,3),X5_5=X(4,4);
  if (X.m()==5)
  {
    return
        X1_1*X2_2*X3_3*X4_4*X5_5 - X1_1*X2_2*X3_3*X4_5*X5_4 - X1_1*X2_2*X3_4*X4_3*X5_5
      + X1_1*X2_2*X3_4*X4_5*X5_3 + X1_1*X2_2*X3_5*X4_3*X5_4 - X1_1*X2_2*X3_5*X4_4*X5_3
      - X1_1*X2_3*X3_2*X4_4*X5_5 + X1_1*X2_3*X3_2*X4_5*X5_4 + X1_1*X2_3*X3_4*X4_2*X5_5
      - X1_1*X2_3*X3_4*X4_5*X5_2 - X1_1*X2_3*X3_5*X4_2*X5_4 + X1_1*X2_3*X3_5*X4_4*X5_2
      + X1_1*X2_4*X3_2*X4_3*X5_5 - X1_1*X2_4*X3_2*X4_5*X5_3 - X1_1*X2_4*X3_3*X4_2*X5_5
      + X1_1*X2_4*X3_3*X4_5*X5_2 + X1_1*X2_4*X3_5*X4_2*X5_3 - X1_1*X2_4*X3_5*X4_3*X5_2
      - X1_1*X2_5*X3_2*X4_3*X5_4 + X1_1*X2_5*X3_2*X4_4*X5_3 + X1_1*X2_5*X3_3*X4_2*X5_4
      - X1_1*X2_5*X3_3*X4_4*X5_2 - X1_1*X2_5*X3_4*X4_2*X5_3 + X1_1*X2_5*X3_4*X4_3*X5_2
      - X1_2*X2_1*X3_3*X4_4*X5_5 + X1_2*X2_1*X3_3*X4_5*X5_4 + X1_2*X2_1*X3_4*X4_3*X5_5
      - X1_2*X2_1*X3_4*X4_5*X5_3 - X1_2*X2_1*X3_5*X4_3*X5_4 + X1_2*X2_1*X3_5*X4_4*X5_3
      + X1_2*X2_3*X3_1*X4_4*X5_5 - X1_2*X2_3*X3_1*X4_5*X5_4 - X1_2*X2_3*X3_4*X4_1*X5_5
      + X1_2*X2_3*X3_4*X4_5*X5_1 + X1_2*X2_3*X3_5*X4_1*X5_4 - X1_2*X2_3*X3_5*X4_4*X5_1
      - X1_2*X2_4*X3_1*X4_3*X5_5 + X1_2*X2_4*X3_1*X4_5*X5_3 + X1_2*X2_4*X3_3*X4_1*X5_5
      - X1_2*X2_4*X3_3*X4_5*X5_1 - X1_2*X2_4*X3_5*X4_1*X5_3 + X1_2*X2_4*X3_5*X4_3*X5_1
      + X1_2*X2_5*X3_1*X4_3*X5_4 - X1_2*X2_5*X3_1*X4_4*X5_3 - X1_2*X2_5*X3_3*X4_1*X5_4
      + X1_2*X2_5*X3_3*X4_4*X5_1 + X1_2*X2_5*X3_4*X4_1*X5_3 - X1_2*X2_5*X3_4*X4_3*X5_1
      + X1_3*X2_1*X3_2*X4_4*X5_5 - X1_3*X2_1*X3_2*X4_5*X5_4 - X1_3*X2_1*X3_4*X4_2*X5_5
      + X1_3*X2_1*X3_4*X4_5*X5_2 + X1_3*X2_1*X3_5*X4_2*X5_4 - X1_3*X2_1*X3_5*X4_4*X5_2
      - X1_3*X2_2*X3_1*X4_4*X5_5 + X1_3*X2_2*X3_1*X4_5*X5_4 + X1_3*X2_2*X3_4*X4_1*X5_5
      - X1_3*X2_2*X3_4*X4_5*X5_1 - X1_3*X2_2*X3_5*X4_1*X5_4 + X1_3*X2_2*X3_5*X4_4*X5_1
      + X1_3*X2_4*X3_1*X4_2*X5_5 - X1_3*X2_4*X3_1*X4_5*X5_2 - X1_3*X2_4*X3_2*X4_1*X5_5
      + X1_3*X2_4*X3_2*X4_5*X5_1 + X1_3*X2_4*X3_5*X4_1*X5_2 - X1_3*X2_4*X3_5*X4_2*X5_1
      - X1_3*X2_5*X3_1*X4_2*X5_4 + X1_3*X2_5*X3_1*X4_4*X5_2 + X1_3*X2_5*X3_2*X4_1*X5_4
      - X1_3*X2_5*X3_2*X4_4*X5_1 - X1_3*X2_5*X3_4*X4_1*X5_2 + X1_3*X2_5*X3_4*X4_2*X5_1
      - X1_4*X2_1*X3_2*X4_3*X5_5 + X1_4*X2_1*X3_2*X4_5*X5_3 + X1_4*X2_1*X3_3*X4_2*X5_5
      - X1_4*X2_1*X3_3*X4_5*X5_2 - X1_4*X2_1*X3_5*X4_2*X5_3 + X1_4*X2_1*X3_5*X4_3*X5_2
      + X1_4*X2_2*X3_1*X4_3*X5_5 - X1_4*X2_2*X3_1*X4_5*X5_3 - X1_4*X2_2*X3_3*X4_1*X5_5
      + X1_4*X2_2*X3_3*X4_5*X5_1 + X1_4*X2_2*X3_5*X4_1*X5_3 - X1_4*X2_2*X3_5*X4_3*X5_1
      - X1_4*X2_3*X3_1*X4_2*X5_5 + X1_4*X2_3*X3_1*X4_5*X5_2 + X1_4*X2_3*X3_2*X4_1*X5_5
      - X1_4*X2_3*X3_2*X4_5*X5_1 - X1_4*X2_3*X3_5*X4_1*X5_2 + X1_4*X2_3*X3_5*X4_2*X5_1
      + X1_4*X2_5*X3_1*X4_2*X5_3 - X1_4*X2_5*X3_1*X4_3*X5_2 - X1_4*X2_5*X3_2*X4_1*X5_3
      + X1_4*X2_5*X3_2*X4_3*X5_1 + X1_4*X2_5*X3_3*X4_1*X5_2 - X1_4*X2_5*X3_3*X4_2*X5_1
      + X1_5*X2_1*X3_2*X4_3*X5_4 - X1_5*X2_1*X3_2*X4_4*X5_3 - X1_5*X2_1*X3_3*X4_2*X5_4
      + X1_5*X2_1*X3_3*X4_4*X5_2 + X1_5*X2_1*X3_4*X4_2*X5_3 - X1_5*X2_1*X3_4*X4_3*X5_2
      - X1_5*X2_2*X3_1*X4_3*X5_4 + X1_5*X2_2*X3_1*X4_4*X5_3 + X1_5*X2_2*X3_3*X4_1*X5_4
      - X1_5*X2_2*X3_3*X4_4*X5_1 - X1_5*X2_2*X3_4*X4_1*X5_3 + X1_5*X2_2*X3_4*X4_3*X5_1
      + X1_5*X2_3*X3_1*X4_2*X5_4 - X1_5*X2_3*X3_1*X4_4*X5_2 - X1_5*X2_3*X3_2*X4_1*X5_4
      + X1_5*X2_3*X3_2*X4_4*X5_1 + X1_5*X2_3*X3_4*X4_1*X5_2 - X1_5*X2_3*X3_4*X4_2*X5_1
      - X1_5*X2_4*X3_1*X4_2*X5_3 + X1_5*X2_4*X3_1*X4_3*X5_2 + X1_5*X2_4*X3_2*X4_1*X5_3
      - X1_5*X2_4*X3_2*X4_3*X5_1 - X1_5*X2_4*X3_3*X4_1*X5_2 + X1_5*X2_4*X3_3*X4_2*X5_1;
  }
  const type X1_6=X(0,5),X2_6=X(1,5),X3_6=X(2,5),X4_6=X(3,5),X5_6=X(4,5);
  const type X6_1=X(5,0),X6_2=X(5,1),X6_3=X(5,2),X6_4=X(5,3),X6_5=X(5,4),X6_6=X(5,5);
  if (X.n()==6)
  {
    return
    X1_1*X2_2*X3_3*X4_4*X5_5*X6_6 - X1_1*X2_2*X3_3*X4_4*X5_6*X6_5 - X1_1*X2_2*X3_3*X4_5*X5_4*X6_6 + X1_1*X2_2*X3_3*X4_5*X5_6*X6_4 + X1_1*X2_2*X3_3*X4_6*X5_4*X6_5 - X1_1*X2_2*X3_3*X4_6*X5_5*X6_4 - X1_1*X2_2*X3_4*X4_3*X5_5*X6_6 + X1_1*X2_2*X3_4*X4_3*X5_6*X6_5 + X1_1*X2_2*X3_4*X4_5*X5_3*X6_6 - X1_1*X2_2*X3_4*X4_5*X5_6*X6_3 - X1_1*X2_2*X3_4*X4_6*X5_3*X6_5 + X1_1*X2_2*X3_4*X4_6*X5_5*X6_3 + X1_1*X2_2*X3_5*X4_3*X5_4*X6_6 - X1_1*X2_2*X3_5*X4_3*X5_6*X6_4 - X1_1*X2_2*X3_5*X4_4*X5_3*X6_6 + X1_1*X2_2*X3_5*X4_4*X5_6*X6_3 + X1_1*X2_2*X3_5*X4_6*X5_3*X6_4 - X1_1*X2_2*X3_5*X4_6*X5_4*X6_3 - X1_1*X2_2*X3_6*X4_3*X5_4*X6_5 + X1_1*X2_2*X3_6*X4_3*X5_5*X6_4 + X1_1*X2_2*X3_6*X4_4*X5_3*X6_5 - X1_1*X2_2*X3_6*X4_4*X5_5*X6_3 - X1_1*X2_2*X3_6*X4_5*X5_3*X6_4 + X1_1*X2_2*X3_6*X4_5*X5_4*X6_3 - X1_1*X2_3*X3_2*X4_4*X5_5*X6_6 + X1_1*X2_3*X3_2*X4_4*X5_6*X6_5 + X1_1*X2_3*X3_2*X4_5*X5_4*X6_6 - X1_1*X2_3*X3_2*X4_5*X5_6*X6_4 - X1_1*X2_3*X3_2*X4_6*X5_4*X6_5 + X1_1*X2_3*X3_2*X4_6*X5_5*X6_4 + X1_1*X2_3*X3_4*X4_2*X5_5*X6_6 - X1_1*X2_3*X3_4*X4_2*X5_6*X6_5 - X1_1*X2_3*X3_4*X4_5*X5_2*X6_6 + X1_1*X2_3*X3_4*X4_5*X5_6*X6_2 + X1_1*X2_3*X3_4*X4_6*X5_2*X6_5 - X1_1*X2_3*X3_4*X4_6*X5_5*X6_2 - X1_1*X2_3*X3_5*X4_2*X5_4*X6_6 + X1_1*X2_3*X3_5*X4_2*X5_6*X6_4 + X1_1*X2_3*X3_5*X4_4*X5_2*X6_6 - X1_1*X2_3*X3_5*X4_4*X5_6*X6_2 - X1_1*X2_3*X3_5*X4_6*X5_2*X6_4 + X1_1*X2_3*X3_5*X4_6*X5_4*X6_2 + X1_1*X2_3*X3_6*X4_2*X5_4*X6_5 - X1_1*X2_3*X3_6*X4_2*X5_5*X6_4 - X1_1*X2_3*X3_6*X4_4*X5_2*X6_5 + X1_1*X2_3*X3_6*X4_4*X5_5*X6_2 + X1_1*X2_3*X3_6*X4_5*X5_2*X6_4 - X1_1*X2_3*X3_6*X4_5*X5_4*X6_2 + X1_1*X2_4*X3_2*X4_3*X5_5*X6_6 - X1_1*X2_4*X3_2*X4_3*X5_6*X6_5 - X1_1*X2_4*X3_2*X4_5*X5_3*X6_6 + X1_1*X2_4*X3_2*X4_5*X5_6*X6_3 + X1_1*X2_4*X3_2*X4_6*X5_3*X6_5 - X1_1*X2_4*X3_2*X4_6*X5_5*X6_3 - X1_1*X2_4*X3_3*X4_2*X5_5*X6_6 + X1_1*X2_4*X3_3*X4_2*X5_6*X6_5 + X1_1*X2_4*X3_3*X4_5*X5_2*X6_6 - X1_1*X2_4*X3_3*X4_5*X5_6*X6_2 - X1_1*X2_4*X3_3*X4_6*X5_2*X6_5 + X1_1*X2_4*X3_3*X4_6*X5_5*X6_2 + X1_1*X2_4*X3_5*X4_2*X5_3*X6_6 - X1_1*X2_4*X3_5*X4_2*X5_6*X6_3 - X1_1*X2_4*X3_5*X4_3*X5_2*X6_6 + X1_1*X2_4*X3_5*X4_3*X5_6*X6_2 + X1_1*X2_4*X3_5*X4_6*X5_2*X6_3 - X1_1*X2_4*X3_5*X4_6*X5_3*X6_2 - X1_1*X2_4*X3_6*X4_2*X5_3*X6_5 + X1_1*X2_4*X3_6*X4_2*X5_5*X6_3 + X1_1*X2_4*X3_6*X4_3*X5_2*X6_5 - X1_1*X2_4*X3_6*X4_3*X5_5*X6_2 - X1_1*X2_4*X3_6*X4_5*X5_2*X6_3 + X1_1*X2_4*X3_6*X4_5*X5_3*X6_2 - X1_1*X2_5*X3_2*X4_3*X5_4*X6_6 + X1_1*X2_5*X3_2*X4_3*X5_6*X6_4 + X1_1*X2_5*X3_2*X4_4*X5_3*X6_6 - X1_1*X2_5*X3_2*X4_4*X5_6*X6_3 - X1_1*X2_5*X3_2*X4_6*X5_3*X6_4 + X1_1*X2_5*X3_2*X4_6*X5_4*X6_3 + X1_1*X2_5*X3_3*X4_2*X5_4*X6_6 - X1_1*X2_5*X3_3*X4_2*X5_6*X6_4 - X1_1*X2_5*X3_3*X4_4*X5_2*X6_6 + X1_1*X2_5*X3_3*X4_4*X5_6*X6_2 + X1_1*X2_5*X3_3*X4_6*X5_2*X6_4 - X1_1*X2_5*X3_3*X4_6*X5_4*X6_2 - X1_1*X2_5*X3_4*X4_2*X5_3*X6_6 + X1_1*X2_5*X3_4*X4_2*X5_6*X6_3 + X1_1*X2_5*X3_4*X4_3*X5_2*X6_6 - X1_1*X2_5*X3_4*X4_3*X5_6*X6_2 - X1_1*X2_5*X3_4*X4_6*X5_2*X6_3 + X1_1*X2_5*X3_4*X4_6*X5_3*X6_2 + X1_1*X2_5*X3_6*X4_2*X5_3*X6_4 - X1_1*X2_5*X3_6*X4_2*X5_4*X6_3 - X1_1*X2_5*X3_6*X4_3*X5_2*X6_4 + X1_1*X2_5*X3_6*X4_3*X5_4*X6_2 + X1_1*X2_5*X3_6*X4_4*X5_2*X6_3 - X1_1*X2_5*X3_6*X4_4*X5_3*X6_2 + X1_1*X2_6*X3_2*X4_3*X5_4*X6_5 - X1_1*X2_6*X3_2*X4_3*X5_5*X6_4 - X1_1*X2_6*X3_2*X4_4*X5_3*X6_5 + X1_1*X2_6*X3_2*X4_4*X5_5*X6_3 + X1_1*X2_6*X3_2*X4_5*X5_3*X6_4 - X1_1*X2_6*X3_2*X4_5*X5_4*X6_3 - X1_1*X2_6*X3_3*X4_2*X5_4*X6_5 + X1_1*X2_6*X3_3*X4_2*X5_5*X6_4 + X1_1*X2_6*X3_3*X4_4*X5_2*X6_5 - X1_1*X2_6*X3_3*X4_4*X5_5*X6_2 - X1_1*X2_6*X3_3*X4_5*X5_2*X6_4 + X1_1*X2_6*X3_3*X4_5*X5_4*X6_2 + X1_1*X2_6*X3_4*X4_2*X5_3*X6_5 - X1_1*X2_6*X3_4*X4_2*X5_5*X6_3 - X1_1*X2_6*X3_4*X4_3*X5_2*X6_5 + X1_1*X2_6*X3_4*X4_3*X5_5*X6_2 + X1_1*X2_6*X3_4*X4_5*X5_2*X6_3 - X1_1*X2_6*X3_4*X4_5*X5_3*X6_2 - X1_1*X2_6*X3_5*X4_2*X5_3*X6_4 + X1_1*X2_6*X3_5*X4_2*X5_4*X6_3 + X1_1*X2_6*X3_5*X4_3*X5_2*X6_4 - X1_1*X2_6*X3_5*X4_3*X5_4*X6_2 - X1_1*X2_6*X3_5*X4_4*X5_2*X6_3 + X1_1*X2_6*X3_5*X4_4*X5_3*X6_2 - X1_2*X2_1*X3_3*X4_4*X5_5*X6_6 + X1_2*X2_1*X3_3*X4_4*X5_6*X6_5 + X1_2*X2_1*X3_3*X4_5*X5_4*X6_6 - X1_2*X2_1*X3_3*X4_5*X5_6*X6_4 - X1_2*X2_1*X3_3*X4_6*X5_4*X6_5 + X1_2*X2_1*X3_3*X4_6*X5_5*X6_4 + X1_2*X2_1*X3_4*X4_3*X5_5*X6_6 - X1_2*X2_1*X3_4*X4_3*X5_6*X6_5 - X1_2*X2_1*X3_4*X4_5*X5_3*X6_6 + X1_2*X2_1*X3_4*X4_5*X5_6*X6_3 + X1_2*X2_1*X3_4*X4_6*X5_3*X6_5 - X1_2*X2_1*X3_4*X4_6*X5_5*X6_3 - X1_2*X2_1*X3_5*X4_3*X5_4*X6_6 + X1_2*X2_1*X3_5*X4_3*X5_6*X6_4 + X1_2*X2_1*X3_5*X4_4*X5_3*X6_6 - X1_2*X2_1*X3_5*X4_4*X5_6*X6_3 - X1_2*X2_1*X3_5*X4_6*X5_3*X6_4 + X1_2*X2_1*X3_5*X4_6*X5_4*X6_3 + X1_2*X2_1*X3_6*X4_3*X5_4*X6_5 - X1_2*X2_1*X3_6*X4_3*X5_5*X6_4 - X1_2*X2_1*X3_6*X4_4*X5_3*X6_5 + X1_2*X2_1*X3_6*X4_4*X5_5*X6_3 + X1_2*X2_1*X3_6*X4_5*X5_3*X6_4 - X1_2*X2_1*X3_6*X4_5*X5_4*X6_3 + X1_2*X2_3*X3_1*X4_4*X5_5*X6_6 - X1_2*X2_3*X3_1*X4_4*X5_6*X6_5 - X1_2*X2_3*X3_1*X4_5*X5_4*X6_6 + X1_2*X2_3*X3_1*X4_5*X5_6*X6_4 + X1_2*X2_3*X3_1*X4_6*X5_4*X6_5 - X1_2*X2_3*X3_1*X4_6*X5_5*X6_4 - X1_2*X2_3*X3_4*X4_1*X5_5*X6_6 + X1_2*X2_3*X3_4*X4_1*X5_6*X6_5 + X1_2*X2_3*X3_4*X4_5*X5_1*X6_6 - X1_2*X2_3*X3_4*X4_5*X5_6*X6_1 - X1_2*X2_3*X3_4*X4_6*X5_1*X6_5 + X1_2*X2_3*X3_4*X4_6*X5_5*X6_1 + X1_2*X2_3*X3_5*X4_1*X5_4*X6_6 - X1_2*X2_3*X3_5*X4_1*X5_6*X6_4 - X1_2*X2_3*X3_5*X4_4*X5_1*X6_6 + X1_2*X2_3*X3_5*X4_4*X5_6*X6_1 + X1_2*X2_3*X3_5*X4_6*X5_1*X6_4 - X1_2*X2_3*X3_5*X4_6*X5_4*X6_1 - X1_2*X2_3*X3_6*X4_1*X5_4*X6_5 + X1_2*X2_3*X3_6*X4_1*X5_5*X6_4 + X1_2*X2_3*X3_6*X4_4*X5_1*X6_5 - X1_2*X2_3*X3_6*X4_4*X5_5*X6_1 - X1_2*X2_3*X3_6*X4_5*X5_1*X6_4 + X1_2*X2_3*X3_6*X4_5*X5_4*X6_1 - X1_2*X2_4*X3_1*X4_3*X5_5*X6_6 + X1_2*X2_4*X3_1*X4_3*X5_6*X6_5 + X1_2*X2_4*X3_1*X4_5*X5_3*X6_6 - X1_2*X2_4*X3_1*X4_5*X5_6*X6_3 - X1_2*X2_4*X3_1*X4_6*X5_3*X6_5 + X1_2*X2_4*X3_1*X4_6*X5_5*X6_3 + X1_2*X2_4*X3_3*X4_1*X5_5*X6_6 - X1_2*X2_4*X3_3*X4_1*X5_6*X6_5 - X1_2*X2_4*X3_3*X4_5*X5_1*X6_6 + X1_2*X2_4*X3_3*X4_5*X5_6*X6_1 + X1_2*X2_4*X3_3*X4_6*X5_1*X6_5 - X1_2*X2_4*X3_3*X4_6*X5_5*X6_1 - X1_2*X2_4*X3_5*X4_1*X5_3*X6_6 + X1_2*X2_4*X3_5*X4_1*X5_6*X6_3 + X1_2*X2_4*X3_5*X4_3*X5_1*X6_6 - X1_2*X2_4*X3_5*X4_3*X5_6*X6_1 - X1_2*X2_4*X3_5*X4_6*X5_1*X6_3 + X1_2*X2_4*X3_5*X4_6*X5_3*X6_1 + X1_2*X2_4*X3_6*X4_1*X5_3*X6_5 - X1_2*X2_4*X3_6*X4_1*X5_5*X6_3 - X1_2*X2_4*X3_6*X4_3*X5_1*X6_5 + X1_2*X2_4*X3_6*X4_3*X5_5*X6_1 + X1_2*X2_4*X3_6*X4_5*X5_1*X6_3 - X1_2*X2_4*X3_6*X4_5*X5_3*X6_1 + X1_2*X2_5*X3_1*X4_3*X5_4*X6_6 - X1_2*X2_5*X3_1*X4_3*X5_6*X6_4 - X1_2*X2_5*X3_1*X4_4*X5_3*X6_6 + X1_2*X2_5*X3_1*X4_4*X5_6*X6_3 + X1_2*X2_5*X3_1*X4_6*X5_3*X6_4 - X1_2*X2_5*X3_1*X4_6*X5_4*X6_3 - X1_2*X2_5*X3_3*X4_1*X5_4*X6_6 + X1_2*X2_5*X3_3*X4_1*X5_6*X6_4 + X1_2*X2_5*X3_3*X4_4*X5_1*X6_6 - X1_2*X2_5*X3_3*X4_4*X5_6*X6_1 - X1_2*X2_5*X3_3*X4_6*X5_1*X6_4 + X1_2*X2_5*X3_3*X4_6*X5_4*X6_1 + X1_2*X2_5*X3_4*X4_1*X5_3*X6_6 - X1_2*X2_5*X3_4*X4_1*X5_6*X6_3 - X1_2*X2_5*X3_4*X4_3*X5_1*X6_6 + X1_2*X2_5*X3_4*X4_3*X5_6*X6_1 + X1_2*X2_5*X3_4*X4_6*X5_1*X6_3 - X1_2*X2_5*X3_4*X4_6*X5_3*X6_1 - X1_2*X2_5*X3_6*X4_1*X5_3*X6_4 + X1_2*X2_5*X3_6*X4_1*X5_4*X6_3 + X1_2*X2_5*X3_6*X4_3*X5_1*X6_4 - X1_2*X2_5*X3_6*X4_3*X5_4*X6_1 - X1_2*X2_5*X3_6*X4_4*X5_1*X6_3 + X1_2*X2_5*X3_6*X4_4*X5_3*X6_1 - X1_2*X2_6*X3_1*X4_3*X5_4*X6_5 + X1_2*X2_6*X3_1*X4_3*X5_5*X6_4 + X1_2*X2_6*X3_1*X4_4*X5_3*X6_5 - X1_2*X2_6*X3_1*X4_4*X5_5*X6_3 - X1_2*X2_6*X3_1*X4_5*X5_3*X6_4 + X1_2*X2_6*X3_1*X4_5*X5_4*X6_3 + X1_2*X2_6*X3_3*X4_1*X5_4*X6_5 - X1_2*X2_6*X3_3*X4_1*X5_5*X6_4 - X1_2*X2_6*X3_3*X4_4*X5_1*X6_5 + X1_2*X2_6*X3_3*X4_4*X5_5*X6_1 + X1_2*X2_6*X3_3*X4_5*X5_1*X6_4 - X1_2*X2_6*X3_3*X4_5*X5_4*X6_1 - X1_2*X2_6*X3_4*X4_1*X5_3*X6_5 + X1_2*X2_6*X3_4*X4_1*X5_5*X6_3 + X1_2*X2_6*X3_4*X4_3*X5_1*X6_5 - X1_2*X2_6*X3_4*X4_3*X5_5*X6_1 - X1_2*X2_6*X3_4*X4_5*X5_1*X6_3 + X1_2*X2_6*X3_4*X4_5*X5_3*X6_1 + X1_2*X2_6*X3_5*X4_1*X5_3*X6_4 - X1_2*X2_6*X3_5*X4_1*X5_4*X6_3 - X1_2*X2_6*X3_5*X4_3*X5_1*X6_4 + X1_2*X2_6*X3_5*X4_3*X5_4*X6_1 + X1_2*X2_6*X3_5*X4_4*X5_1*X6_3 - X1_2*X2_6*X3_5*X4_4*X5_3*X6_1 + X1_3*X2_1*X3_2*X4_4*X5_5*X6_6 - X1_3*X2_1*X3_2*X4_4*X5_6*X6_5 - X1_3*X2_1*X3_2*X4_5*X5_4*X6_6 + X1_3*X2_1*X3_2*X4_5*X5_6*X6_4 + X1_3*X2_1*X3_2*X4_6*X5_4*X6_5 - X1_3*X2_1*X3_2*X4_6*X5_5*X6_4 - X1_3*X2_1*X3_4*X4_2*X5_5*X6_6 + X1_3*X2_1*X3_4*X4_2*X5_6*X6_5 + X1_3*X2_1*X3_4*X4_5*X5_2*X6_6 - X1_3*X2_1*X3_4*X4_5*X5_6*X6_2 - X1_3*X2_1*X3_4*X4_6*X5_2*X6_5 + X1_3*X2_1*X3_4*X4_6*X5_5*X6_2 + X1_3*X2_1*X3_5*X4_2*X5_4*X6_6 - X1_3*X2_1*X3_5*X4_2*X5_6*X6_4 - X1_3*X2_1*X3_5*X4_4*X5_2*X6_6 + X1_3*X2_1*X3_5*X4_4*X5_6*X6_2 + X1_3*X2_1*X3_5*X4_6*X5_2*X6_4 - X1_3*X2_1*X3_5*X4_6*X5_4*X6_2 - X1_3*X2_1*X3_6*X4_2*X5_4*X6_5 + X1_3*X2_1*X3_6*X4_2*X5_5*X6_4 + X1_3*X2_1*X3_6*X4_4*X5_2*X6_5 - X1_3*X2_1*X3_6*X4_4*X5_5*X6_2 - X1_3*X2_1*X3_6*X4_5*X5_2*X6_4 + X1_3*X2_1*X3_6*X4_5*X5_4*X6_2 - X1_3*X2_2*X3_1*X4_4*X5_5*X6_6 + X1_3*X2_2*X3_1*X4_4*X5_6*X6_5 + X1_3*X2_2*X3_1*X4_5*X5_4*X6_6 - X1_3*X2_2*X3_1*X4_5*X5_6*X6_4 - X1_3*X2_2*X3_1*X4_6*X5_4*X6_5 + X1_3*X2_2*X3_1*X4_6*X5_5*X6_4 + X1_3*X2_2*X3_4*X4_1*X5_5*X6_6 - X1_3*X2_2*X3_4*X4_1*X5_6*X6_5 - X1_3*X2_2*X3_4*X4_5*X5_1*X6_6 + X1_3*X2_2*X3_4*X4_5*X5_6*X6_1 + X1_3*X2_2*X3_4*X4_6*X5_1*X6_5 - X1_3*X2_2*X3_4*X4_6*X5_5*X6_1 - X1_3*X2_2*X3_5*X4_1*X5_4*X6_6 + X1_3*X2_2*X3_5*X4_1*X5_6*X6_4 + X1_3*X2_2*X3_5*X4_4*X5_1*X6_6 - X1_3*X2_2*X3_5*X4_4*X5_6*X6_1 - X1_3*X2_2*X3_5*X4_6*X5_1*X6_4 + X1_3*X2_2*X3_5*X4_6*X5_4*X6_1 + X1_3*X2_2*X3_6*X4_1*X5_4*X6_5 - X1_3*X2_2*X3_6*X4_1*X5_5*X6_4 - X1_3*X2_2*X3_6*X4_4*X5_1*X6_5 + X1_3*X2_2*X3_6*X4_4*X5_5*X6_1 + X1_3*X2_2*X3_6*X4_5*X5_1*X6_4 - X1_3*X2_2*X3_6*X4_5*X5_4*X6_1 + X1_3*X2_4*X3_1*X4_2*X5_5*X6_6 - X1_3*X2_4*X3_1*X4_2*X5_6*X6_5 - X1_3*X2_4*X3_1*X4_5*X5_2*X6_6 + X1_3*X2_4*X3_1*X4_5*X5_6*X6_2 + X1_3*X2_4*X3_1*X4_6*X5_2*X6_5 - X1_3*X2_4*X3_1*X4_6*X5_5*X6_2 - X1_3*X2_4*X3_2*X4_1*X5_5*X6_6 + X1_3*X2_4*X3_2*X4_1*X5_6*X6_5 + X1_3*X2_4*X3_2*X4_5*X5_1*X6_6 - X1_3*X2_4*X3_2*X4_5*X5_6*X6_1 - X1_3*X2_4*X3_2*X4_6*X5_1*X6_5 + X1_3*X2_4*X3_2*X4_6*X5_5*X6_1 + X1_3*X2_4*X3_5*X4_1*X5_2*X6_6 - X1_3*X2_4*X3_5*X4_1*X5_6*X6_2 - X1_3*X2_4*X3_5*X4_2*X5_1*X6_6 + X1_3*X2_4*X3_5*X4_2*X5_6*X6_1 + X1_3*X2_4*X3_5*X4_6*X5_1*X6_2 - X1_3*X2_4*X3_5*X4_6*X5_2*X6_1 - X1_3*X2_4*X3_6*X4_1*X5_2*X6_5 + X1_3*X2_4*X3_6*X4_1*X5_5*X6_2 + X1_3*X2_4*X3_6*X4_2*X5_1*X6_5 - X1_3*X2_4*X3_6*X4_2*X5_5*X6_1 - X1_3*X2_4*X3_6*X4_5*X5_1*X6_2 + X1_3*X2_4*X3_6*X4_5*X5_2*X6_1 - X1_3*X2_5*X3_1*X4_2*X5_4*X6_6 + X1_3*X2_5*X3_1*X4_2*X5_6*X6_4 + X1_3*X2_5*X3_1*X4_4*X5_2*X6_6 - X1_3*X2_5*X3_1*X4_4*X5_6*X6_2 - X1_3*X2_5*X3_1*X4_6*X5_2*X6_4 + X1_3*X2_5*X3_1*X4_6*X5_4*X6_2 + X1_3*X2_5*X3_2*X4_1*X5_4*X6_6 - X1_3*X2_5*X3_2*X4_1*X5_6*X6_4 - X1_3*X2_5*X3_2*X4_4*X5_1*X6_6 + X1_3*X2_5*X3_2*X4_4*X5_6*X6_1 + X1_3*X2_5*X3_2*X4_6*X5_1*X6_4 - X1_3*X2_5*X3_2*X4_6*X5_4*X6_1 - X1_3*X2_5*X3_4*X4_1*X5_2*X6_6 + X1_3*X2_5*X3_4*X4_1*X5_6*X6_2 + X1_3*X2_5*X3_4*X4_2*X5_1*X6_6 - X1_3*X2_5*X3_4*X4_2*X5_6*X6_1 - X1_3*X2_5*X3_4*X4_6*X5_1*X6_2 + X1_3*X2_5*X3_4*X4_6*X5_2*X6_1 + X1_3*X2_5*X3_6*X4_1*X5_2*X6_4 - X1_3*X2_5*X3_6*X4_1*X5_4*X6_2 - X1_3*X2_5*X3_6*X4_2*X5_1*X6_4 + X1_3*X2_5*X3_6*X4_2*X5_4*X6_1 + X1_3*X2_5*X3_6*X4_4*X5_1*X6_2 - X1_3*X2_5*X3_6*X4_4*X5_2*X6_1 + X1_3*X2_6*X3_1*X4_2*X5_4*X6_5 - X1_3*X2_6*X3_1*X4_2*X5_5*X6_4 - X1_3*X2_6*X3_1*X4_4*X5_2*X6_5 + X1_3*X2_6*X3_1*X4_4*X5_5*X6_2 + X1_3*X2_6*X3_1*X4_5*X5_2*X6_4 - X1_3*X2_6*X3_1*X4_5*X5_4*X6_2 - X1_3*X2_6*X3_2*X4_1*X5_4*X6_5 + X1_3*X2_6*X3_2*X4_1*X5_5*X6_4 + X1_3*X2_6*X3_2*X4_4*X5_1*X6_5 - X1_3*X2_6*X3_2*X4_4*X5_5*X6_1 - X1_3*X2_6*X3_2*X4_5*X5_1*X6_4 + X1_3*X2_6*X3_2*X4_5*X5_4*X6_1 + X1_3*X2_6*X3_4*X4_1*X5_2*X6_5 - X1_3*X2_6*X3_4*X4_1*X5_5*X6_2 - X1_3*X2_6*X3_4*X4_2*X5_1*X6_5 + X1_3*X2_6*X3_4*X4_2*X5_5*X6_1 + X1_3*X2_6*X3_4*X4_5*X5_1*X6_2 - X1_3*X2_6*X3_4*X4_5*X5_2*X6_1 - X1_3*X2_6*X3_5*X4_1*X5_2*X6_4 + X1_3*X2_6*X3_5*X4_1*X5_4*X6_2 + X1_3*X2_6*X3_5*X4_2*X5_1*X6_4 - X1_3*X2_6*X3_5*X4_2*X5_4*X6_1 - X1_3*X2_6*X3_5*X4_4*X5_1*X6_2 + X1_3*X2_6*X3_5*X4_4*X5_2*X6_1 - X1_4*X2_1*X3_2*X4_3*X5_5*X6_6 + X1_4*X2_1*X3_2*X4_3*X5_6*X6_5 + X1_4*X2_1*X3_2*X4_5*X5_3*X6_6 - X1_4*X2_1*X3_2*X4_5*X5_6*X6_3 - X1_4*X2_1*X3_2*X4_6*X5_3*X6_5 + X1_4*X2_1*X3_2*X4_6*X5_5*X6_3 + X1_4*X2_1*X3_3*X4_2*X5_5*X6_6 - X1_4*X2_1*X3_3*X4_2*X5_6*X6_5 - X1_4*X2_1*X3_3*X4_5*X5_2*X6_6 + X1_4*X2_1*X3_3*X4_5*X5_6*X6_2 + X1_4*X2_1*X3_3*X4_6*X5_2*X6_5 - X1_4*X2_1*X3_3*X4_6*X5_5*X6_2 - X1_4*X2_1*X3_5*X4_2*X5_3*X6_6 + X1_4*X2_1*X3_5*X4_2*X5_6*X6_3 + X1_4*X2_1*X3_5*X4_3*X5_2*X6_6 - X1_4*X2_1*X3_5*X4_3*X5_6*X6_2 - X1_4*X2_1*X3_5*X4_6*X5_2*X6_3 + X1_4*X2_1*X3_5*X4_6*X5_3*X6_2 + X1_4*X2_1*X3_6*X4_2*X5_3*X6_5 - X1_4*X2_1*X3_6*X4_2*X5_5*X6_3 - X1_4*X2_1*X3_6*X4_3*X5_2*X6_5 + X1_4*X2_1*X3_6*X4_3*X5_5*X6_2 + X1_4*X2_1*X3_6*X4_5*X5_2*X6_3 - X1_4*X2_1*X3_6*X4_5*X5_3*X6_2 + X1_4*X2_2*X3_1*X4_3*X5_5*X6_6 - X1_4*X2_2*X3_1*X4_3*X5_6*X6_5 - X1_4*X2_2*X3_1*X4_5*X5_3*X6_6 + X1_4*X2_2*X3_1*X4_5*X5_6*X6_3 + X1_4*X2_2*X3_1*X4_6*X5_3*X6_5 - X1_4*X2_2*X3_1*X4_6*X5_5*X6_3 - X1_4*X2_2*X3_3*X4_1*X5_5*X6_6 + X1_4*X2_2*X3_3*X4_1*X5_6*X6_5 + X1_4*X2_2*X3_3*X4_5*X5_1*X6_6 - X1_4*X2_2*X3_3*X4_5*X5_6*X6_1 - X1_4*X2_2*X3_3*X4_6*X5_1*X6_5 + X1_4*X2_2*X3_3*X4_6*X5_5*X6_1 + X1_4*X2_2*X3_5*X4_1*X5_3*X6_6 - X1_4*X2_2*X3_5*X4_1*X5_6*X6_3 - X1_4*X2_2*X3_5*X4_3*X5_1*X6_6 + X1_4*X2_2*X3_5*X4_3*X5_6*X6_1 + X1_4*X2_2*X3_5*X4_6*X5_1*X6_3 - X1_4*X2_2*X3_5*X4_6*X5_3*X6_1 - X1_4*X2_2*X3_6*X4_1*X5_3*X6_5 + X1_4*X2_2*X3_6*X4_1*X5_5*X6_3 + X1_4*X2_2*X3_6*X4_3*X5_1*X6_5 - X1_4*X2_2*X3_6*X4_3*X5_5*X6_1 - X1_4*X2_2*X3_6*X4_5*X5_1*X6_3 + X1_4*X2_2*X3_6*X4_5*X5_3*X6_1 - X1_4*X2_3*X3_1*X4_2*X5_5*X6_6 + X1_4*X2_3*X3_1*X4_2*X5_6*X6_5 + X1_4*X2_3*X3_1*X4_5*X5_2*X6_6 - X1_4*X2_3*X3_1*X4_5*X5_6*X6_2 - X1_4*X2_3*X3_1*X4_6*X5_2*X6_5 + X1_4*X2_3*X3_1*X4_6*X5_5*X6_2 + X1_4*X2_3*X3_2*X4_1*X5_5*X6_6 - X1_4*X2_3*X3_2*X4_1*X5_6*X6_5 - X1_4*X2_3*X3_2*X4_5*X5_1*X6_6 + X1_4*X2_3*X3_2*X4_5*X5_6*X6_1 + X1_4*X2_3*X3_2*X4_6*X5_1*X6_5 - X1_4*X2_3*X3_2*X4_6*X5_5*X6_1 - X1_4*X2_3*X3_5*X4_1*X5_2*X6_6 + X1_4*X2_3*X3_5*X4_1*X5_6*X6_2 + X1_4*X2_3*X3_5*X4_2*X5_1*X6_6 - X1_4*X2_3*X3_5*X4_2*X5_6*X6_1 - X1_4*X2_3*X3_5*X4_6*X5_1*X6_2 + X1_4*X2_3*X3_5*X4_6*X5_2*X6_1 + X1_4*X2_3*X3_6*X4_1*X5_2*X6_5 - X1_4*X2_3*X3_6*X4_1*X5_5*X6_2 - X1_4*X2_3*X3_6*X4_2*X5_1*X6_5 + X1_4*X2_3*X3_6*X4_2*X5_5*X6_1 + X1_4*X2_3*X3_6*X4_5*X5_1*X6_2 - X1_4*X2_3*X3_6*X4_5*X5_2*X6_1 + X1_4*X2_5*X3_1*X4_2*X5_3*X6_6 - X1_4*X2_5*X3_1*X4_2*X5_6*X6_3 - X1_4*X2_5*X3_1*X4_3*X5_2*X6_6 + X1_4*X2_5*X3_1*X4_3*X5_6*X6_2 + X1_4*X2_5*X3_1*X4_6*X5_2*X6_3 - X1_4*X2_5*X3_1*X4_6*X5_3*X6_2 - X1_4*X2_5*X3_2*X4_1*X5_3*X6_6 + X1_4*X2_5*X3_2*X4_1*X5_6*X6_3 + X1_4*X2_5*X3_2*X4_3*X5_1*X6_6 - X1_4*X2_5*X3_2*X4_3*X5_6*X6_1 - X1_4*X2_5*X3_2*X4_6*X5_1*X6_3 + X1_4*X2_5*X3_2*X4_6*X5_3*X6_1 + X1_4*X2_5*X3_3*X4_1*X5_2*X6_6 - X1_4*X2_5*X3_3*X4_1*X5_6*X6_2 - X1_4*X2_5*X3_3*X4_2*X5_1*X6_6 + X1_4*X2_5*X3_3*X4_2*X5_6*X6_1 + X1_4*X2_5*X3_3*X4_6*X5_1*X6_2 - X1_4*X2_5*X3_3*X4_6*X5_2*X6_1 - X1_4*X2_5*X3_6*X4_1*X5_2*X6_3 + X1_4*X2_5*X3_6*X4_1*X5_3*X6_2 + X1_4*X2_5*X3_6*X4_2*X5_1*X6_3 - X1_4*X2_5*X3_6*X4_2*X5_3*X6_1 - X1_4*X2_5*X3_6*X4_3*X5_1*X6_2 + X1_4*X2_5*X3_6*X4_3*X5_2*X6_1 - X1_4*X2_6*X3_1*X4_2*X5_3*X6_5 + X1_4*X2_6*X3_1*X4_2*X5_5*X6_3 + X1_4*X2_6*X3_1*X4_3*X5_2*X6_5 - X1_4*X2_6*X3_1*X4_3*X5_5*X6_2 - X1_4*X2_6*X3_1*X4_5*X5_2*X6_3 + X1_4*X2_6*X3_1*X4_5*X5_3*X6_2 + X1_4*X2_6*X3_2*X4_1*X5_3*X6_5 - X1_4*X2_6*X3_2*X4_1*X5_5*X6_3 - X1_4*X2_6*X3_2*X4_3*X5_1*X6_5 + X1_4*X2_6*X3_2*X4_3*X5_5*X6_1 + X1_4*X2_6*X3_2*X4_5*X5_1*X6_3 - X1_4*X2_6*X3_2*X4_5*X5_3*X6_1 - X1_4*X2_6*X3_3*X4_1*X5_2*X6_5 + X1_4*X2_6*X3_3*X4_1*X5_5*X6_2 + X1_4*X2_6*X3_3*X4_2*X5_1*X6_5 - X1_4*X2_6*X3_3*X4_2*X5_5*X6_1 - X1_4*X2_6*X3_3*X4_5*X5_1*X6_2 + X1_4*X2_6*X3_3*X4_5*X5_2*X6_1 + X1_4*X2_6*X3_5*X4_1*X5_2*X6_3 - X1_4*X2_6*X3_5*X4_1*X5_3*X6_2 - X1_4*X2_6*X3_5*X4_2*X5_1*X6_3 + X1_4*X2_6*X3_5*X4_2*X5_3*X6_1 + X1_4*X2_6*X3_5*X4_3*X5_1*X6_2 - X1_4*X2_6*X3_5*X4_3*X5_2*X6_1 + X1_5*X2_1*X3_2*X4_3*X5_4*X6_6 - X1_5*X2_1*X3_2*X4_3*X5_6*X6_4 - X1_5*X2_1*X3_2*X4_4*X5_3*X6_6 + X1_5*X2_1*X3_2*X4_4*X5_6*X6_3 + X1_5*X2_1*X3_2*X4_6*X5_3*X6_4 - X1_5*X2_1*X3_2*X4_6*X5_4*X6_3 - X1_5*X2_1*X3_3*X4_2*X5_4*X6_6 + X1_5*X2_1*X3_3*X4_2*X5_6*X6_4 + X1_5*X2_1*X3_3*X4_4*X5_2*X6_6 - X1_5*X2_1*X3_3*X4_4*X5_6*X6_2 - X1_5*X2_1*X3_3*X4_6*X5_2*X6_4 + X1_5*X2_1*X3_3*X4_6*X5_4*X6_2 + X1_5*X2_1*X3_4*X4_2*X5_3*X6_6 - X1_5*X2_1*X3_4*X4_2*X5_6*X6_3 - X1_5*X2_1*X3_4*X4_3*X5_2*X6_6 + X1_5*X2_1*X3_4*X4_3*X5_6*X6_2 + X1_5*X2_1*X3_4*X4_6*X5_2*X6_3 - X1_5*X2_1*X3_4*X4_6*X5_3*X6_2 - X1_5*X2_1*X3_6*X4_2*X5_3*X6_4 + X1_5*X2_1*X3_6*X4_2*X5_4*X6_3 + X1_5*X2_1*X3_6*X4_3*X5_2*X6_4 - X1_5*X2_1*X3_6*X4_3*X5_4*X6_2 - X1_5*X2_1*X3_6*X4_4*X5_2*X6_3 + X1_5*X2_1*X3_6*X4_4*X5_3*X6_2 - X1_5*X2_2*X3_1*X4_3*X5_4*X6_6 + X1_5*X2_2*X3_1*X4_3*X5_6*X6_4 + X1_5*X2_2*X3_1*X4_4*X5_3*X6_6 - X1_5*X2_2*X3_1*X4_4*X5_6*X6_3 - X1_5*X2_2*X3_1*X4_6*X5_3*X6_4 + X1_5*X2_2*X3_1*X4_6*X5_4*X6_3 + X1_5*X2_2*X3_3*X4_1*X5_4*X6_6 - X1_5*X2_2*X3_3*X4_1*X5_6*X6_4 - X1_5*X2_2*X3_3*X4_4*X5_1*X6_6 + X1_5*X2_2*X3_3*X4_4*X5_6*X6_1 + X1_5*X2_2*X3_3*X4_6*X5_1*X6_4 - X1_5*X2_2*X3_3*X4_6*X5_4*X6_1 - X1_5*X2_2*X3_4*X4_1*X5_3*X6_6 + X1_5*X2_2*X3_4*X4_1*X5_6*X6_3 + X1_5*X2_2*X3_4*X4_3*X5_1*X6_6 - X1_5*X2_2*X3_4*X4_3*X5_6*X6_1 - X1_5*X2_2*X3_4*X4_6*X5_1*X6_3 + X1_5*X2_2*X3_4*X4_6*X5_3*X6_1 + X1_5*X2_2*X3_6*X4_1*X5_3*X6_4 - X1_5*X2_2*X3_6*X4_1*X5_4*X6_3 - X1_5*X2_2*X3_6*X4_3*X5_1*X6_4 + X1_5*X2_2*X3_6*X4_3*X5_4*X6_1 + X1_5*X2_2*X3_6*X4_4*X5_1*X6_3 - X1_5*X2_2*X3_6*X4_4*X5_3*X6_1 + X1_5*X2_3*X3_1*X4_2*X5_4*X6_6 - X1_5*X2_3*X3_1*X4_2*X5_6*X6_4 - X1_5*X2_3*X3_1*X4_4*X5_2*X6_6 + X1_5*X2_3*X3_1*X4_4*X5_6*X6_2 + X1_5*X2_3*X3_1*X4_6*X5_2*X6_4 - X1_5*X2_3*X3_1*X4_6*X5_4*X6_2 - X1_5*X2_3*X3_2*X4_1*X5_4*X6_6 + X1_5*X2_3*X3_2*X4_1*X5_6*X6_4 + X1_5*X2_3*X3_2*X4_4*X5_1*X6_6 - X1_5*X2_3*X3_2*X4_4*X5_6*X6_1 - X1_5*X2_3*X3_2*X4_6*X5_1*X6_4 + X1_5*X2_3*X3_2*X4_6*X5_4*X6_1 + X1_5*X2_3*X3_4*X4_1*X5_2*X6_6 - X1_5*X2_3*X3_4*X4_1*X5_6*X6_2 - X1_5*X2_3*X3_4*X4_2*X5_1*X6_6 + X1_5*X2_3*X3_4*X4_2*X5_6*X6_1 + X1_5*X2_3*X3_4*X4_6*X5_1*X6_2 - X1_5*X2_3*X3_4*X4_6*X5_2*X6_1 - X1_5*X2_3*X3_6*X4_1*X5_2*X6_4 + X1_5*X2_3*X3_6*X4_1*X5_4*X6_2 + X1_5*X2_3*X3_6*X4_2*X5_1*X6_4 - X1_5*X2_3*X3_6*X4_2*X5_4*X6_1 - X1_5*X2_3*X3_6*X4_4*X5_1*X6_2 + X1_5*X2_3*X3_6*X4_4*X5_2*X6_1 - X1_5*X2_4*X3_1*X4_2*X5_3*X6_6 + X1_5*X2_4*X3_1*X4_2*X5_6*X6_3 + X1_5*X2_4*X3_1*X4_3*X5_2*X6_6 - X1_5*X2_4*X3_1*X4_3*X5_6*X6_2 - X1_5*X2_4*X3_1*X4_6*X5_2*X6_3 + X1_5*X2_4*X3_1*X4_6*X5_3*X6_2 + X1_5*X2_4*X3_2*X4_1*X5_3*X6_6 - X1_5*X2_4*X3_2*X4_1*X5_6*X6_3 - X1_5*X2_4*X3_2*X4_3*X5_1*X6_6 + X1_5*X2_4*X3_2*X4_3*X5_6*X6_1 + X1_5*X2_4*X3_2*X4_6*X5_1*X6_3 - X1_5*X2_4*X3_2*X4_6*X5_3*X6_1 - X1_5*X2_4*X3_3*X4_1*X5_2*X6_6 + X1_5*X2_4*X3_3*X4_1*X5_6*X6_2 + X1_5*X2_4*X3_3*X4_2*X5_1*X6_6 - X1_5*X2_4*X3_3*X4_2*X5_6*X6_1 - X1_5*X2_4*X3_3*X4_6*X5_1*X6_2 + X1_5*X2_4*X3_3*X4_6*X5_2*X6_1 + X1_5*X2_4*X3_6*X4_1*X5_2*X6_3 - X1_5*X2_4*X3_6*X4_1*X5_3*X6_2 - X1_5*X2_4*X3_6*X4_2*X5_1*X6_3 + X1_5*X2_4*X3_6*X4_2*X5_3*X6_1 + X1_5*X2_4*X3_6*X4_3*X5_1*X6_2 - X1_5*X2_4*X3_6*X4_3*X5_2*X6_1 + X1_5*X2_6*X3_1*X4_2*X5_3*X6_4 - X1_5*X2_6*X3_1*X4_2*X5_4*X6_3 - X1_5*X2_6*X3_1*X4_3*X5_2*X6_4 + X1_5*X2_6*X3_1*X4_3*X5_4*X6_2 + X1_5*X2_6*X3_1*X4_4*X5_2*X6_3 - X1_5*X2_6*X3_1*X4_4*X5_3*X6_2 - X1_5*X2_6*X3_2*X4_1*X5_3*X6_4 + X1_5*X2_6*X3_2*X4_1*X5_4*X6_3 + X1_5*X2_6*X3_2*X4_3*X5_1*X6_4 - X1_5*X2_6*X3_2*X4_3*X5_4*X6_1 - X1_5*X2_6*X3_2*X4_4*X5_1*X6_3 + X1_5*X2_6*X3_2*X4_4*X5_3*X6_1 + X1_5*X2_6*X3_3*X4_1*X5_2*X6_4 - X1_5*X2_6*X3_3*X4_1*X5_4*X6_2 - X1_5*X2_6*X3_3*X4_2*X5_1*X6_4 + X1_5*X2_6*X3_3*X4_2*X5_4*X6_1 + X1_5*X2_6*X3_3*X4_4*X5_1*X6_2 - X1_5*X2_6*X3_3*X4_4*X5_2*X6_1 - X1_5*X2_6*X3_4*X4_1*X5_2*X6_3 + X1_5*X2_6*X3_4*X4_1*X5_3*X6_2 + X1_5*X2_6*X3_4*X4_2*X5_1*X6_3 - X1_5*X2_6*X3_4*X4_2*X5_3*X6_1 - X1_5*X2_6*X3_4*X4_3*X5_1*X6_2 + X1_5*X2_6*X3_4*X4_3*X5_2*X6_1 - X1_6*X2_1*X3_2*X4_3*X5_4*X6_5 + X1_6*X2_1*X3_2*X4_3*X5_5*X6_4 + X1_6*X2_1*X3_2*X4_4*X5_3*X6_5 - X1_6*X2_1*X3_2*X4_4*X5_5*X6_3 - X1_6*X2_1*X3_2*X4_5*X5_3*X6_4 + X1_6*X2_1*X3_2*X4_5*X5_4*X6_3 + X1_6*X2_1*X3_3*X4_2*X5_4*X6_5 - X1_6*X2_1*X3_3*X4_2*X5_5*X6_4 - X1_6*X2_1*X3_3*X4_4*X5_2*X6_5 + X1_6*X2_1*X3_3*X4_4*X5_5*X6_2 + X1_6*X2_1*X3_3*X4_5*X5_2*X6_4 - X1_6*X2_1*X3_3*X4_5*X5_4*X6_2 - X1_6*X2_1*X3_4*X4_2*X5_3*X6_5 + X1_6*X2_1*X3_4*X4_2*X5_5*X6_3 + X1_6*X2_1*X3_4*X4_3*X5_2*X6_5 - X1_6*X2_1*X3_4*X4_3*X5_5*X6_2 - X1_6*X2_1*X3_4*X4_5*X5_2*X6_3 + X1_6*X2_1*X3_4*X4_5*X5_3*X6_2 + X1_6*X2_1*X3_5*X4_2*X5_3*X6_4 - X1_6*X2_1*X3_5*X4_2*X5_4*X6_3 - X1_6*X2_1*X3_5*X4_3*X5_2*X6_4 + X1_6*X2_1*X3_5*X4_3*X5_4*X6_2 + X1_6*X2_1*X3_5*X4_4*X5_2*X6_3 - X1_6*X2_1*X3_5*X4_4*X5_3*X6_2 + X1_6*X2_2*X3_1*X4_3*X5_4*X6_5 - X1_6*X2_2*X3_1*X4_3*X5_5*X6_4 - X1_6*X2_2*X3_1*X4_4*X5_3*X6_5 + X1_6*X2_2*X3_1*X4_4*X5_5*X6_3 + X1_6*X2_2*X3_1*X4_5*X5_3*X6_4 - X1_6*X2_2*X3_1*X4_5*X5_4*X6_3 - X1_6*X2_2*X3_3*X4_1*X5_4*X6_5 + X1_6*X2_2*X3_3*X4_1*X5_5*X6_4 + X1_6*X2_2*X3_3*X4_4*X5_1*X6_5 - X1_6*X2_2*X3_3*X4_4*X5_5*X6_1 - X1_6*X2_2*X3_3*X4_5*X5_1*X6_4 + X1_6*X2_2*X3_3*X4_5*X5_4*X6_1 + X1_6*X2_2*X3_4*X4_1*X5_3*X6_5 - X1_6*X2_2*X3_4*X4_1*X5_5*X6_3 - X1_6*X2_2*X3_4*X4_3*X5_1*X6_5 + X1_6*X2_2*X3_4*X4_3*X5_5*X6_1 + X1_6*X2_2*X3_4*X4_5*X5_1*X6_3 - X1_6*X2_2*X3_4*X4_5*X5_3*X6_1 - X1_6*X2_2*X3_5*X4_1*X5_3*X6_4 + X1_6*X2_2*X3_5*X4_1*X5_4*X6_3 + X1_6*X2_2*X3_5*X4_3*X5_1*X6_4 - X1_6*X2_2*X3_5*X4_3*X5_4*X6_1 - X1_6*X2_2*X3_5*X4_4*X5_1*X6_3 + X1_6*X2_2*X3_5*X4_4*X5_3*X6_1 - X1_6*X2_3*X3_1*X4_2*X5_4*X6_5 + X1_6*X2_3*X3_1*X4_2*X5_5*X6_4 + X1_6*X2_3*X3_1*X4_4*X5_2*X6_5 - X1_6*X2_3*X3_1*X4_4*X5_5*X6_2 - X1_6*X2_3*X3_1*X4_5*X5_2*X6_4 + X1_6*X2_3*X3_1*X4_5*X5_4*X6_2 + X1_6*X2_3*X3_2*X4_1*X5_4*X6_5 - X1_6*X2_3*X3_2*X4_1*X5_5*X6_4 - X1_6*X2_3*X3_2*X4_4*X5_1*X6_5 + X1_6*X2_3*X3_2*X4_4*X5_5*X6_1 + X1_6*X2_3*X3_2*X4_5*X5_1*X6_4 - X1_6*X2_3*X3_2*X4_5*X5_4*X6_1 - X1_6*X2_3*X3_4*X4_1*X5_2*X6_5 + X1_6*X2_3*X3_4*X4_1*X5_5*X6_2 + X1_6*X2_3*X3_4*X4_2*X5_1*X6_5 - X1_6*X2_3*X3_4*X4_2*X5_5*X6_1 - X1_6*X2_3*X3_4*X4_5*X5_1*X6_2 + X1_6*X2_3*X3_4*X4_5*X5_2*X6_1 + X1_6*X2_3*X3_5*X4_1*X5_2*X6_4 - X1_6*X2_3*X3_5*X4_1*X5_4*X6_2 - X1_6*X2_3*X3_5*X4_2*X5_1*X6_4 + X1_6*X2_3*X3_5*X4_2*X5_4*X6_1 + X1_6*X2_3*X3_5*X4_4*X5_1*X6_2 - X1_6*X2_3*X3_5*X4_4*X5_2*X6_1 + X1_6*X2_4*X3_1*X4_2*X5_3*X6_5 - X1_6*X2_4*X3_1*X4_2*X5_5*X6_3 - X1_6*X2_4*X3_1*X4_3*X5_2*X6_5 + X1_6*X2_4*X3_1*X4_3*X5_5*X6_2 + X1_6*X2_4*X3_1*X4_5*X5_2*X6_3 - X1_6*X2_4*X3_1*X4_5*X5_3*X6_2 - X1_6*X2_4*X3_2*X4_1*X5_3*X6_5 + X1_6*X2_4*X3_2*X4_1*X5_5*X6_3 + X1_6*X2_4*X3_2*X4_3*X5_1*X6_5 - X1_6*X2_4*X3_2*X4_3*X5_5*X6_1 - X1_6*X2_4*X3_2*X4_5*X5_1*X6_3 + X1_6*X2_4*X3_2*X4_5*X5_3*X6_1 + X1_6*X2_4*X3_3*X4_1*X5_2*X6_5 - X1_6*X2_4*X3_3*X4_1*X5_5*X6_2 - X1_6*X2_4*X3_3*X4_2*X5_1*X6_5 + X1_6*X2_4*X3_3*X4_2*X5_5*X6_1 + X1_6*X2_4*X3_3*X4_5*X5_1*X6_2 - X1_6*X2_4*X3_3*X4_5*X5_2*X6_1 - X1_6*X2_4*X3_5*X4_1*X5_2*X6_3 + X1_6*X2_4*X3_5*X4_1*X5_3*X6_2 + X1_6*X2_4*X3_5*X4_2*X5_1*X6_3 - X1_6*X2_4*X3_5*X4_2*X5_3*X6_1 - X1_6*X2_4*X3_5*X4_3*X5_1*X6_2 + X1_6*X2_4*X3_5*X4_3*X5_2*X6_1 - X1_6*X2_5*X3_1*X4_2*X5_3*X6_4 + X1_6*X2_5*X3_1*X4_2*X5_4*X6_3 + X1_6*X2_5*X3_1*X4_3*X5_2*X6_4 - X1_6*X2_5*X3_1*X4_3*X5_4*X6_2 - X1_6*X2_5*X3_1*X4_4*X5_2*X6_3 + X1_6*X2_5*X3_1*X4_4*X5_3*X6_2 + X1_6*X2_5*X3_2*X4_1*X5_3*X6_4 - X1_6*X2_5*X3_2*X4_1*X5_4*X6_3 - X1_6*X2_5*X3_2*X4_3*X5_1*X6_4 + X1_6*X2_5*X3_2*X4_3*X5_4*X6_1 + X1_6*X2_5*X3_2*X4_4*X5_1*X6_3 - X1_6*X2_5*X3_2*X4_4*X5_3*X6_1 - X1_6*X2_5*X3_3*X4_1*X5_2*X6_4 + X1_6*X2_5*X3_3*X4_1*X5_4*X6_2 + X1_6*X2_5*X3_3*X4_2*X5_1*X6_4 - X1_6*X2_5*X3_3*X4_2*X5_4*X6_1 - X1_6*X2_5*X3_3*X4_4*X5_1*X6_2 + X1_6*X2_5*X3_3*X4_4*X5_2*X6_1 + X1_6*X2_5*X3_4*X4_1*X5_2*X6_3 - X1_6*X2_5*X3_4*X4_1*X5_3*X6_2 - X1_6*X2_5*X3_4*X4_2*X5_1*X6_3 + X1_6*X2_5*X3_4*X4_2*X5_3*X6_1 + X1_6*X2_5*X3_4*X4_3*X5_1*X6_2 - X1_6*X2_5*X3_4*X4_3*X5_2*X6_1;

  }
  // compute the eigenvalues
  vecd<type> lambda( X.m() );
  eign( X , lambda );
  type d = 1;
  for (index_t i = 0; i < X.m(); i++)
    d *= lambda[i];
  return d;
  printf("unsupported n = %d\n",int(X.m()));
  avro_assert(false);
}

template<typename type>
type
det(const symd<type>& X)
{
  avro_assert(X.m()==X.n());
  if (X.m() == 1) return X(0,0);
  if (X.m() == 2) return X(1,1)*X(0,0)-X(0,1)*X(1,0);
  if (X.m() == 3) return X(0,0)*(X(2,2)*X(1,1)-X(2,1)*X(1,2))
                    -X(1,0)*(X(2,2)*X(0,1)-X(2,1)*X(0,2))
                    +X(2,0)*(X(1,2)*X(0,1)-X(1,1)*X(0,2));
  const type X1_1 = X(0,0),X1_2=X(0,1),X1_3=X(0,2),X1_4=X(0,3);
  const type X2_1 = X(1,0),X2_2=X(1,1),X2_3=X(1,2),X2_4=X(1,3);
  const type X3_1 = X(2,0),X3_2=X(2,1),X3_3=X(2,2),X3_4=X(2,3);
  const type X4_1 = X(3,0),X4_2=X(3,1),X4_3=X(3,2),X4_4=X(3,3);
  if (X.m() == 4)
  {   return
        X1_1*X2_2*X3_3*X4_4 - X1_1*X2_2*X3_4*X4_3 - X1_1*X2_3*X3_2*X4_4
      + X1_1*X2_3*X3_4*X4_2 + X1_1*X2_4*X3_2*X4_3 - X1_1*X2_4*X3_3*X4_2
      - X1_2*X2_1*X3_3*X4_4 + X1_2*X2_1*X3_4*X4_3 + X1_2*X2_3*X3_1*X4_4
      - X1_2*X2_3*X3_4*X4_1 - X1_2*X2_4*X3_1*X4_3 + X1_2*X2_4*X3_3*X4_1
      + X1_3*X2_1*X3_2*X4_4 - X1_3*X2_1*X3_4*X4_2 - X1_3*X2_2*X3_1*X4_4
      + X1_3*X2_2*X3_4*X4_1 + X1_3*X2_4*X3_1*X4_2 - X1_3*X2_4*X3_2*X4_1
      - X1_4*X2_1*X3_2*X4_3 + X1_4*X2_1*X3_3*X4_2 + X1_4*X2_2*X3_1*X4_3
      - X1_4*X2_2*X3_3*X4_1 - X1_4*X2_3*X3_1*X4_2 + X1_4*X2_3*X3_2*X4_1;
  }
  printf("unsupported n = %d\n",int(X.m()));
  avro_assert(false);
}

template<typename type>
symd<type>
expm( const symd<type>& m ) {
	std::pair< vecd<type>,matd<type> > decomp = m.eig();
	for (index_t k = 0; k < m.n(); k++)
		decomp.first(k) = ::exp( decomp.first(k) );
	return symd<type>(decomp);
}

template<typename type>
symd<type>
logm( const symd<type>& m ) {
  std::pair< vecd<type>,matd<type> > decomp = m.eig();
	for (index_t k = 0; k < m.n(); k++)
		decomp.first(k) = ::log( decomp.first(k) );
	return symd<type>(decomp);
}

template<typename type>
symd<type>
powm( const symd<type>& m , real_t p ) {
  std::pair< vecd<type>,matd<type> > decomp = m.eig();
	for (index_t k = 0; k < m.n(); k++)
		decomp.first(k) = ::pow( decomp.first(k) , p );
	return symd<type>(decomp);
}

template<typename type>
symd<type>
sqrtm( const symd<type>& m ) {
  std::pair< vecd<type>,matd<type> > decomp = m.eig();
	for (index_t k = 0; k < m.n(); k++)
		decomp.first(k) = ::sqrt( decomp.first(k) );
	return symd<type>(decomp);
}

template<typename type>
type
quadratic_form( const symd<type>& M , const vecd<real_t>& e ) {
	if (M.n() == 2) {
		type mu = M(0,0)*e(0) + M(0,1)*e(1);
		type mv = M(1,0)*e(0) + M(1,1)*e(1);
		return mu*e(0) + mv*e(1);
	}
	else if (M.n() == 3) {
		type mu = M(0,0)*e(0) + M(0,1)*e(1) + M(0,2)*e(2);
		type mv = M(1,0)*e(0) + M(1,1)*e(1) + M(1,2)*e(2);
		type mw = M(2,0)*e(0) + M(2,1)*e(1) + M(2,2)*e(2);
		return mu*e(0) + mv*e(1) + mw*e(2);
	}
	else if (M.n() == 4) {
		type mu = M(0,0)*e(0) + M(0,1)*e(1) + M(0,2)*e(2) + M(0,3)*e(3);
		type mv = M(1,0)*e(0) + M(1,1)*e(1) + M(1,2)*e(2) + M(1,3)*e(3);
		type mw = M(2,0)*e(0) + M(2,1)*e(1) + M(2,2)*e(2) + M(2,3)*e(3);
		type mt = M(3,0)*e(0) + M(3,1)*e(1) + M(3,2)*e(2) + M(3,3)*e(3);
		return mu*e(0) + mv*e(1) + mw*e(2) + mt*e(3);
	}
	else
		avro_implement;
	return -1.;
}

template void solveLUP( const matd<real_t>& , const vecd<real_t>& , vecd<real_t>& );
template void inverseLUP( const matd<real_t>& , matd<real_t>& );

#define INSTANTIATE_TRANSPOSE(T) template matd<T> transpose( const matd<T>& );
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

#define INSTANTIATE_INV(T) \
  template matd<T> inverse( const matd<T>& ); \
  template symd<T> inverse( const symd<T>& );
INSTANTIATE_INV( real_t )
template matd<float> inverse( const matd<float>& );
#undef INSTANTIATE_INV

#define INSTANTIATE_EIG(T) \
  template void eig( const symd<T>& m , vecd<T>& L , matd<T>& Q ); \
  template std::pair< vecd<T> , matd<T> > eig( const symd<T>& m );
INSTANTIATE_EIG( real_t )
INSTANTIATE_EIG( SurrealS<1> );
#undef INSTANTIATE_EIG

// compiling takes really long with surreals
template<>
SurrealD
det(const matd<SurrealD>& X) {
  avro_implement;
}

#define INST_DETERMINANT(X) \
  template<> SurrealS<X> det(const matd<SurrealS<X>>&) \
    { avro_implement; }

#if USE_SURREAL
INST_DETERMINANT(1)
INST_DETERMINANT(2)
INST_DETERMINANT(3)
INST_DETERMINANT(4)
INST_DETERMINANT(6)
INST_DETERMINANT(10)
#endif

template real_t det(const matd<real_t>& X);
template real_t det(const symd<real_t>& X);
template dual   det(const matd<dual>& X);

#define INSTANTIATE_SYMD_FUNC(T) \
  template symd<T> logm( const symd<T>& ); \
  template symd<T> expm( const symd<T>& ); \
  template symd<T> sqrtm( const symd<T>& ); \
  template symd<T> powm( const symd<T>& , real_t ); \
  template symd<T> interp( const std::vector<real_t>& , const std::vector<symd<T>>& ); \
  template T quadratic_form( const symd<T>& , const vecd<real_t>& );
INSTANTIATE_SYMD_FUNC( real_t )
INSTANTIATE_SYMD_FUNC( SurrealS<1> )
INSTANTIATE_SYMD_FUNC( SurrealS<3> )
INSTANTIATE_SYMD_FUNC( SurrealS<6> )
INSTANTIATE_SYMD_FUNC( SurrealS<10> )
#undef INSTANTIATE_SYMD_FUNC

} // numerics

} // avro
