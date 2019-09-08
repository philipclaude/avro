#ifndef URSA_LIB_NUMERICS_COORDINATE_H_
#define URSA_LIB_NUMERICS_COORDINATE_H_

#include "common/types.h"

#include <vector>

namespace ursa
{

namespace numerics
{

class Coordinate : public std::vector<real>
{
public:
  using std::vector<real>::data;
  using std::vector<real>::operator[];

  Coordinate( const coord_t dim ) :
    std::vector<real>(dim)
  {}

  Coordinate( real* x , const coord_t dim ) :
    std::vector<real>(x,x+dim)
  {}

  coord_t dim() const { return size(); }

private:
  using std::vector<real>::size;

};

} // numerics

} // ursa

#endif
