#ifndef URSA_LIB_LIBRARY_OBJ_H_
#define URSA_LIB_LIBRARY_OBJ_H_

#include "mesh/topology.h"
#include "mesh/points.h"

#include <string>

namespace ursa
{

namespace library
{

class objFile : public Topology<Simplex<Lagrange>>
{
public:
  objFile( const std::string& filename );

  void read();

private:
  Points points_;
  std::string filename_;
};

} // library

} // ursa

#endif
