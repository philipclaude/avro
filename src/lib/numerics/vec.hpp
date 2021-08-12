#include "numerics/vec.h"

namespace avro
{

template<typename R, typename S>
vecd< typename result_of<R,S>::type >
operator+ ( const vecd<R>& x , const vecd<S>& y ) {
  typedef typename result_of<R,S>::type T;
  avro_assert( x.m() == x.m() );
  vecd<T> z( x.m() );
  for (index_t i = 0; i < z.m(); i++)
    z(i) = x(i) + y(i);
  return z;
}

template<typename R, typename S>
vecd< typename result_of<R,S>::type >
operator- ( const vecd<R>& x , const vecd<S>& y ) {
  typedef typename result_of<R,S>::type T;
  avro_assert( x.m() == x.m() );
  vecd<T> z( x.m() );
  for (index_t i = 0; i < z.m(); i++)
    z(i) = x(i) - y(i);
  return z;
}


#define INSTANTIATE_VECADD(R,S,T) \
template<index_t M> \
vecs<M,T> \
operator+ ( const vecs<M,R>& u , const vecs<M,S>& v ) { \
  vecs<M,T> w; \
  for (index_t i = 0; i < M; i++) \
    w(i) = u(i) + v(i); \
  return w; \
}

#define INSTANTIATE_VECSUB(R,S,T) \
template<index_t M> \
vecs<M,T> \
operator- ( const vecs<M,R>& u , const vecs<M,S>& v ) { \
  vecs<M,T> w; \
  for (index_t i = 0; i < M; i++) \
    w(i) = u(i) - v(i); \
  return w; \
}

#define INSTANTIATE_VECINC(S,T) \
template<index_t M> \
vecs<M,T>& \
operator+= ( vecs<M,S>& u , const vecs<M,T>& v ) { \
  for (index_t i = 0; i < M; i++) \
    u(i) += v(i); \
  return u; \
}

#define INSTANTIATE_VECDEC(S,T) \
template<index_t M> \
vecs<M,S>& \
operator-= ( vecs<M,S>& u , const vecs<M,T>& v ) { \
  for (index_t i = 0; i < M; i++) \
    u(i) -= v(i); \
  return u; \
}

#define INSTANTIATE_VECVECMUL(R,S,T) \
template<index_t M> \
vecs<M,T> operator* (const vecs<M,R>& u, const vecs<M,S>& v) { \
  vecs<M,T> w; \
  for (index_t i = 0; i < M ; i++) \
    w(i) = u(i)*v(i); \
  return w; \
}

#define INSTANTIATE_VECSCAMUL_R(R,S,T) \
template<index_t M> \
vecs<M,T> operator* (const vecs<M,R>& u, const S& a) { \
  vecs<M,T> v; \
  for (index_t i = 0; i < M ; i++) \
    v(i) = a*u(i); \
  return v; \
}

#define INSTANTIATE_VECSCAMUL_L(R,S,T) \
template<index_t M> \
vecs<M,T> operator* (const R& a, const vecs<M,S>& u) { \
  vecs<M,T> v; \
  for (index_t i = 0; i < M ; i++) \
    v(i) = a*u(i); \
  return v; \
}

#define INSTANTIATE_DOT(R,S,T) \
template<index_t M> \
T \
dot( const vecs<M,R>& u , const vecs<M,S>& v ) { \
  T result = 0; \
  for (index_t i = 0; i < M; i++) \
     result += u(i)*v(i); \
  return result; \
}

template<index_t M, typename T>
T norm( const vecs<M,T>& u ) {
  return sqrt( dot(u,u) );
}

namespace numerics
{
#define INSTANTIATE_NORMALIZE(T) \
template<index_t M> \
void normalize( vecs<M,T>& u ) { \
  T n = norm(u); \
  for (index_t i = 0; i < M; i++) \
    u(i) /= n; \
}

} // numerics

#define COMMA ,

INSTANTIATE_VECADD( real_t , real_t , real_t )
INSTANTIATE_VECADD( float , float , float )
INSTANTIATE_VECSUB( real_t , real_t , real_t )
INSTANTIATE_VECSUB( float , float , float )
INSTANTIATE_VECADD( vecs<2 COMMA real_t> , vecs<2 COMMA real_t> , vecs<2 COMMA real_t> )
INSTANTIATE_VECADD( vecs<3 COMMA real_t> , vecs<3 COMMA real_t> , vecs<3 COMMA real_t> )

INSTANTIATE_VECINC( real_t , real_t )
INSTANTIATE_VECDEC( real_t , real_t )

INSTANTIATE_VECINC( vecs<2 COMMA real_t > , vecs<2 COMMA real_t> )
INSTANTIATE_VECINC( vecs<3 COMMA real_t > , vecs<3 COMMA real_t> )
INSTANTIATE_VECINC( SurrealS<1> , SurrealS<1> )

INSTANTIATE_VECINC( SurrealS<2> , SurrealS<2> )
INSTANTIATE_VECDEC( SurrealS<2> , SurrealS<2> )
INSTANTIATE_VECDEC( SurrealS<2> , real_t )

INSTANTIATE_VECINC( SurrealS<3> , SurrealS<3> )
INSTANTIATE_VECDEC( SurrealS<3> , SurrealS<3> )
INSTANTIATE_VECDEC( SurrealS<3> , real_t )

INSTANTIATE_VECVECMUL( real_t , vecs< 2 COMMA SurrealS<2> > , vecs< 2 COMMA SurrealS<2> > )
INSTANTIATE_VECVECMUL( real_t , vecs< 3 COMMA SurrealS<3> > , vecs< 3 COMMA SurrealS<3> > )

INSTANTIATE_DOT( real_t , real_t , real_t )
INSTANTIATE_DOT( SurrealS<1> , SurrealS<1> , SurrealS<1> )

INSTANTIATE_DOT( real_t , SurrealS<1> , SurrealS<1> )
INSTANTIATE_DOT( real_t , SurrealS<2> , SurrealS<2> )
INSTANTIATE_DOT( real_t , SurrealS<3> , SurrealS<3> )

INSTANTIATE_DOT( real_t , vecs<2 COMMA SurrealS<2> > , vecs<2 COMMA SurrealS<2> > )
INSTANTIATE_DOT( real_t , vecs<3 COMMA SurrealS<3> > , vecs<3 COMMA SurrealS<3> > )

INSTANTIATE_VECSCAMUL_R( real_t , real_t , real_t )
INSTANTIATE_VECSCAMUL_R( float , float , float )
INSTANTIATE_VECSCAMUL_R( vecs<2 COMMA real_t> , real_t , vecs<2 COMMA real_t> )
INSTANTIATE_VECSCAMUL_R( vecs<3 COMMA real_t> , real_t , vecs<3 COMMA real_t> )

INSTANTIATE_VECSCAMUL_L( real_t , real_t , real_t )
INSTANTIATE_VECSCAMUL_L( float , float , float )
INSTANTIATE_VECSCAMUL_L( real_t , SurrealS<1> , SurrealS<1> )
INSTANTIATE_VECSCAMUL_L( real_t , SurrealS<2> , SurrealS<2> )
INSTANTIATE_VECSCAMUL_L( real_t , SurrealS<3> , SurrealS<3> )

namespace numerics
{
INSTANTIATE_NORMALIZE( real_t )
}

#undef COMMA

}
