#ifndef URSA_LIB_LIBRARY_MESHB_H_
#define URSA_LIB_LIBRARY_MESHB_H_

#include "mesh/topology.h"

#include <string>

namespace ursa
{

namespace library
{

class meshbFile : public Topology<Simplex<Lagrange>>
{

public:
  meshbFile( const std::string& filename );

private:
  std::string filename_;
  Points points_;
};

} // library

} // ursa

#endif
