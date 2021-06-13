#include "numerics/sym.h"
#include "numerics/mat.hpp"
#include "numerics/vec.h"

namespace avro
{

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

template vecd<real_t> operator*( const matd<real_t>& , const vecd<real_t>& );

#define INSTANTIATE_MATD(T) template class matd<T>;
INSTANTIATE_MATD(real_t)
#undef INSTANTIATE_MATD

#define INSTANTIATE_MUL(R,S) template symd< typename result_of<R,S>::type > operator*( const symd<R>& , const symd<S>& );
INSTANTIATE_MUL( real_t , real_t )
INSTANTIATE_MUL( SurrealS<1> , SurrealS<1> )
INSTANTIATE_MUL( SurrealS<3> , SurrealS<3> )
INSTANTIATE_MUL( SurrealS<6> , SurrealS<6> )
INSTANTIATE_MUL( SurrealS<10> , SurrealS<10> )
#undef INSTANTIATE_MUL

#define INSTANTIATE_MULD(R,S) template matd< typename result_of<R,S>::type > operator*( const matd<R>& , const matd<S>& );
INSTANTIATE_MULD( real_t , real_t )
INSTANTIATE_MULD( SurrealS<1> , SurrealS<1> )
INSTANTIATE_MULD( SurrealS<3> , SurrealS<3> )
INSTANTIATE_MULD( SurrealS<6> , SurrealS<6> )
INSTANTIATE_MULD( SurrealS<10> , SurrealS<10> )
#undef INSTANTIATE_MULD

#define INSTANTIATE_OPD(R,S) \
  template matd< typename result_of<R,S>::type > operator+( const matd<R>& , const matd<S>& ); \
  template matd< typename result_of<R,S>::type > operator-( const matd<R>& , const matd<S>& );
INSTANTIATE_OPD( real_t , real_t )
INSTANTIATE_OPD( SurrealS<1> , SurrealS<1> )
INSTANTIATE_OPD( SurrealS<3> , SurrealS<3> )
INSTANTIATE_OPD( SurrealS<6> , SurrealS<6> )
INSTANTIATE_OPD( SurrealS<10> , SurrealS<10> )
#undef INSTANTIATE_OPD


} // avro
