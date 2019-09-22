#ifndef URSA_LIB_GRAPHICS_PLOTTER_H_
#define URSA_LIB_GRAPHICS_PLOTTER_H_

#include "graphics/client.h"
#include "graphics/server.h"

#include <vector>

namespace ursa
{

namespace graphics
{

template<typename Framework>
class Plotter
{



private:
  Server<Framework> server_;
  std::vector< Client<Framework> > client_;

};


} // graphics

} // ursa

#endif
