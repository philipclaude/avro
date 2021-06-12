#include "common/error.h"
#include "numerics/linear_algebra.h"
#include "numerics/mat.h"
#include "numerics/vec.h"

namespace avro
{

template<typename T>
void
matd<T>::set_row( index_t i , const vecd<T>& row ) {
  avro_assert( row.m() == n_ );
  for (index_t j = 0; j < row.m(); j++)
    (*this)(i,j) = row(j);
}

template<typename T>
void
matd<T>::get_row( index_t i , vecd<T>& row ) const {
  avro_assert( row.m() == n_ );
  for (index_t j = 0; j < row.m(); j++)
    row(j) = (*this)(i,j);
}

template<typename R, typename S>
matd< typename result_of<R,S>::type >
operator* (const matd<R>& A, const matd<S>& B) {
  typedef typename result_of<R,S>::type T;
  avro_assert_msg( A.n() == B.m() , "bad matrix sizes" );
  matd<T> C( A.m() , B.n() );
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

template<typename R, typename S>
matd< typename result_of<R,S>::type >
operator+ (const matd<R>& A, const matd<S>& B) {
  typedef typename result_of<R,S>::type T;
  avro_assert_msg( A.m() == B.m() , "bad matrix sizes" );
  avro_assert_msg( A.n() == B.n() , "bad matrix sizes" );
  matd<T> C( A.m() , A.n() );
  for (index_t i = 0; i < C.m(); i++) {
    for (index_t j = 0; j < C.n(); j++) {
      C(i,j) = A(i,j) - B(i,j);
    }
  }
  return C;
}

template<typename R, typename S>
matd< typename result_of<R,S>::type >
operator- (const matd<R>& A, const matd<S>& B) {
  typedef typename result_of<R,S>::type T;
  avro_assert_msg( A.m() == B.m() , "bad matrix sizes" );
  avro_assert_msg( A.n() == B.n() , "bad matrix sizes" );
  matd<T> C( A.m() , A.n() );
  for (index_t i = 0; i < C.m(); i++) {
    for (index_t j = 0; j < C.n(); j++) {
      C(i,j) = A(i,j) - B(i,j);
    }
  }
  return C;
}

template<typename R, typename S>
symd< typename result_of<R,S>::type >
operator* (const symd<R>& A, const symd<S>& B) {
  typedef typename result_of<R,S>::type T;
  avro_assert_msg( A.n() == B.m() , "bad matrix sizes" );
  symd<T> C( A.m() , B.n() );
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

#define INSTANTIATE_MATADD(R,S,T) \
template<index_t M,index_t N> \
mats<M,N,T> \
operator+ ( const mats<M,N,R>& A , const mats<M,N,S>& B ) { \
  mats<M,N,T> C; \
  for (index_t i = 0; i < M; i++) \
  for (index_t j = 0; j < N; j++) \
    C(i,j) = A(i,j) + B(i,j); \
  return C; \
}

#define INSTANTIATE_MATSUB(R,S,T) \
template<index_t M,index_t N> \
mats<M,N,T> \
operator- ( const mats<M,N,R>& A , const mats<M,N,S>& B ) { \
  mats<M,N,T> C; \
  for (index_t i = 0; i < M; i++) \
  for (index_t j = 0; j < N; j++) \
    C(i,j) = A(i,j) + B(i,j); \
  return C; \
}

#define INSTANTIATE_MATMUL(R,S,T) \
template<index_t MA, index_t NA , index_t MB,index_t NB> \
mats<MA,NB,T> operator* (const mats<MA,NA,R>& A, const mats<MB,NB,S>& B) { \
  static_assert( NA == MB , "bad matrix sizes" ); \
  mats<MA,NB,T> C; \
  for (index_t i = 0; i < MA; i++) { \
    for (index_t j = 0; j < NB; j++) { \
      T sum = 0; \
      for (index_t k = 0; k < NA; k++) \
        sum += A(i,k)*B(k,j); \
      C(i,j) = sum; \
    } \
  } \
  return C; \
}

#define INSTANTIATE_MATSCAMUL_L(R,S,T) \
template<index_t M, index_t N> \
mats<M,N,T> operator* (const mats<M,N,R>& A, const S& b) { \
  mats<M,N,T> C; \
  for (index_t i = 0; i < M; i++) \
  for (index_t j = 0; j < N; j++) \
    C(i,j) = A(i,j)*b; \
  return C; \
}

#define INSTANTIATE_MATSCAMUL_R(R,S,T) \
template<index_t M, index_t N> \
mats<M,N,T> operator* ( const R& b , const mats<M,N,S>& A ) { \
  mats<M,N,T> C; \
  for (index_t i = 0; i < M; i++) \
  for (index_t j = 0; j < N; j++) \
    C(i,j) = b*A(i,j); \
  return C; \
}

#define INSTANTIATE_MATVECMUL( R , S ) \
template<index_t M, index_t N> \
vecs< M , typename result_of<R,S>::type > operator*( const mats<M,N,R>& A , const vecs<N,S>& b ) { \
  typedef typename result_of<R,S>::type T; \
  vecs<M,T> c; \
  for (index_t i = 0; i < M; i++) { \
    c(i) = 0; \
    for (index_t j = 0; j < N; j++) \
      c(i) += A(i,j) * b(j); \
  } \
  return c; \
}

namespace numerics
{

template<index_t M, index_t N, typename T>
mats<N,M,T>
transpose( const mats<M,N,T>& A ) {
  mats<N,M,T> At;
  for (index_t i = 0; i < M; i++)
  for (index_t j = 0; j < N; j++)
    At(j,i) = A(i,j);
  return At;
}

template<typename T>
T
det( const mats<2,2,T>& A ) {
  return A(0,0)*A(1,1) - A(0,1)*A(1,0);
}

template<typename T>
T
det( const mats<3,3,T>& A ) {
  return A(0,0)*(A(2,2)*A(1,1)-A(2,1)*A(1,2))
                -A(1,0)*(A(2,2)*A(0,1)-A(2,1)*A(0,2))
                +A(2,0)*(A(1,2)*A(0,1)-A(1,1)*A(0,2));
}

template<typename T>
mats<2,2,T>
inverse( const mats<2,2,T>& A ) {
  mats<2,2,T> Ainv;
  const T id = 1.0/det(A);
  Ainv(0,0) =  A(1,1)*id;
  Ainv(0,1) = -A(0,1)*id;
  Ainv(1,0) = -A(1,0)*id;
  Ainv(1,1) =  A(0,0)*id;
  return Ainv;
}

template<index_t N, typename T>
mats<N,N,T>
inverse( const mats<N,N,T>& A ) {
  matd<T> a(N,N);
  for (index_t i = 0; i < N; i++)
  for (index_t j = 0; j < N; j++)
    a(i,j) = A(i,j);
  matd<T> ai = numerics::inverse(a);
  mats<N,N,T> Ai;
  for (index_t i = 0; i < N; i++)
  for (index_t j = 0; j < N; j++)
    Ai(i,j) = ai(i,j);
  return Ai;
}

} // numerics

#define COMMA ,

INSTANTIATE_MATADD( real_t , real_t , real_t )
INSTANTIATE_MATADD( real_t , SurrealS<4> , SurrealS<4> )
INSTANTIATE_MATADD( SurrealS<4> , real_t , SurrealS<4> )
INSTANTIATE_MATSUB( real_t , real_t , real_t )
INSTANTIATE_MATSUB( real_t , SurrealS<4> , SurrealS<4> )
INSTANTIATE_MATSUB( SurrealS<4> , real_t , SurrealS<4> )

INSTANTIATE_MATMUL( real_t , real_t , real_t )
INSTANTIATE_MATMUL( float , float , float )
INSTANTIATE_MATMUL( real_t , SurrealS<4> , SurrealS<4> )
INSTANTIATE_MATMUL( SurrealS<4> , real_t , SurrealS<4> )

INSTANTIATE_MATVECMUL( float , float )

INSTANTIATE_MATSCAMUL_R( real_t , SurrealS<4> , SurrealS<4> )
//INSTANTIATE_MATSCAMUL_R( real_t , real_t , SurrealS<4> )
INSTANTIATE_MATSCAMUL_L( real_t , SurrealS<4> , SurrealS<4> )
//INSTANTIATE_MATSCAMUL_L( real_t , real_t , SurrealS<4> )
INSTANTIATE_MATSCAMUL_L( real_t , real_t , real_t )
INSTANTIATE_MATSCAMUL_R( real_t , real_t , real_t )

#undef COMMA

} // avro
