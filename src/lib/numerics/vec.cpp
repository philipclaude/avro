#include "vec.hpp"

namespace avro
{

#define INSTANTIATE_VECD(R,S) \
  template vecd< typename result_of<R,S>::type > operator+ (const vecd<R>& , const vecd<S>& ); \
  template vecd< typename result_of<R,S>::type > operator- (const vecd<R>& , const vecd<S>& );

INSTANTIATE_VECD( real_t , real_t )

#undef INSTANTIATE_VECD

namespace numerics
{

template<typename T>
vecs<3,T>
cross( const vecs<3,T>& u , const vecs<3,T>& v ) {
  vecs<3,T> w;
  w(0) =    u(1)*v(2) - u(2)*v(1);
  w(1) = -( u(0)*v(2) - u(2)*v(0) );
  w(2) =    u(0)*v(1) - u(1)*v(0);
  return w;
}

template vecs<3,real_t> cross( const vecs<3,real_t>& , const vecs<3,real_t>& );

} // numerics

} // avro
