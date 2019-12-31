#ifndef avro_LIB_MESH_TRIANGULATION_H_
#define avro_LIB_MESH_TRIANGULATION_H_

#include "common/types.h"

#include "master/simplex.h"
#include "mesh/topology.h"

namespace avro
{

template<typename type>
class Triangulation : public Topology<Simplex>
{
public:
  Triangulation( const Topology<type>& topology );

  void extract();

private:
  const Topology<type>& topology_;
  Points points_;
  std::vector<index_t> parents_;
};

} // avro

#endif
