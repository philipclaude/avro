#include "voronoi/new/facets.h"

namespace avro
{

namespace voronoi
{

bool
operator==( const Bisector& bx , const Bisector& by ) {
  if (bx.p0 != by.p0) return false;
  if (bx.p1 != by.p1) return false;
  return true;
}

// needed to create a map of elements
bool
operator<( const Bisector& f , const Bisector& g ) {
  if (f.p0 < g.p0) return true;
  if (f.p0 == g.p0) {
    if (f.p1 < g.p1) return true;
  }
  return false;
}

} // voronoi

} // avro
