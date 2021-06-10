#include "mat.hpp"
#include "vec.h"

namespace avro
{

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

    for (int i = n-1; i >= 0; i--) {
      for (int k = i+1; k < n; k++)
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

template<typename T>
matd<T>
operator* (const matd<T>& A, const matd<T>& B) {
  avro_assert_msg( A.n() == B.m() , "bad matrix sizes" );
  matd<T> C(A.m(),B.n());
  for (index_t i = 0; i < A.m(); i++) {
    for (index_t j = 0; j < B.n(); j++) {
      T sum = 0;
      for (index_t k = 0; k < A.n(); k++)
        sum += A(i,k)*B(k,j);
      C(i,j) = sum;
    }
  }
  return C;
}

template<typename T>
vecd<T>
operator* (const matd<T>& A, const vecd<T>& x) {
  avro_assert_msg( A.n() == x.m() , "bad sizes" );
  vecd<T> b(A.m());
  for (index_t i = 0; i < A.m(); i++) {
    T sum = 0;
    for (index_t k = 0; k < A.n(); k++)
      sum += A(i,k)*x(k);
    b(i) = sum;
  }
  return b;
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

template void solveLUP( const matd<real_t>& , const vecd<real_t>& , vecd<real_t>& );
template void inverseLUP( const matd<real_t>& , matd<real_t>& );
template vecd<real_t> operator*( const matd<real_t>& , const vecd<real_t>& );
template matd<real_t> operator*( const matd<real_t>& , const matd<real_t>& );

} // avro
