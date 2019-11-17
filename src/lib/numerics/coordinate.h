#ifndef LUNA_LIB_NUMERICS_COORDINATE_H_
#define LUNA_LIB_NUMERICS_COORDINATE_H_

#include "common/types.h"

#include <vector>

namespace luna
{

namespace numerics
{

class Coordinate : public std::vector<real_t>
{
public:
  using std::vector<real_t>::data;
  using std::vector<real_t>::operator[];

  Coordinate( const coord_t dim ) :
    std::vector<real_t>(dim)
  {}

  Coordinate(real_t* x , const coord_t dim ) :
    std::vector<real_t>(x,x+dim)
  {}

  coord_t dim() const { return size(); }

private:
  using std::vector<real_t>::size;

};

} // numerics

} // luna

#endif
