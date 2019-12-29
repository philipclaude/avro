#ifndef luma_LIB_LIBRARY_EGADS_H_
#define luma_LIB_LIBRARY_EGADS_H_

#include "geometry/egads/body.h"

namespace luma
{

namespace EGADS
{

class Cube : public Body
{
public:
  Cube( const Context* context , coord_t dim );
  Cube( const Context* context , const std::vector<real_t>& lens , const real_t* x0=nullptr );

private:
  ego object_;
};

} // EGADS

} // luma

#endif
