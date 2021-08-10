#ifndef AVRO_LIB_VORONOI_NEW_FACETS_H_
#define AVRO_LIB_VORONOI_NEW_FACETS_H_

#include "avro_types.h"

namespace avro
{

namespace voronoi
{

struct Bisector {
  Bisector( index_t q0 , index_t q1 ) {
    if (q0 < q1) {
      p0 = q0;
      p1 = q1;
    }
    else {
      p0 = q1;
      p1 = q0;
    }
  }
  int p0;
  int p1;
};

// needed to create a set/map of bisectors
bool operator==( const Bisector& bx , const Bisector& by );
bool operator< ( const Bisector& f  , const Bisector& g );

} // voronoi

} // avro

#endif
