#ifndef LUNA_LIB_GEOMETRY_MODEL_H_
#define LUNA_LIB_GEOMETRY_MODEL_H_

#include "common/error.h"

#include "geometry/body.h"

namespace luna
{

class Model
{
public:
  Model( coord_t number ) :
    number_(number)
  {}

protected:
  coord_t number_;
};

} // luna

#endif
