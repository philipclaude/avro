#ifndef luma_LIB_LIBRARY_OBJ_H_
#define luma_LIB_LIBRARY_OBJ_H_

#include "mesh/topology.h"
#include "mesh/points.h"

#include <string>

namespace luma
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

} // luma

#endif
