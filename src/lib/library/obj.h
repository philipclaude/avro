#ifndef LUNA_LIB_LIBRARY_OBJ_H_
#define LUNA_LIB_LIBRARY_OBJ_H_

#include "mesh/topology.h"
#include "mesh/points.h"

#include <string>

namespace luna
{

namespace library
{

class objFile : public Topology<Simplex>
{
public:
  objFile( const std::string& filename );

  void read();

private:
  Points points_;
  std::string filename_;
};

} // library

} // luna

#endif
