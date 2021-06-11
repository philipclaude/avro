#include "vec.hpp"

namespace avro
{

#define INSTANTIATE_VECD(R,S) \
  template vecd< typename result_of<R,S>::type > operator+ (const vecd<R>& , const vecd<S>& ); \
  template vecd< typename result_of<R,S>::type > operator- (const vecd<R>& , const vecd<S>& );

INSTANTIATE_VECD( real_t , real_t )

#undef INSTANTIATE_VECD

}