#ifndef avro_LIB_MESH_TRIANGULATION_H_
#define avro_LIB_MESH_TRIANGULATION_H_

#include "common/types.h"

#include "master/simplex.h"
#include "mesh/topology.h"

namespace avro
{

class TriangulationBase : public Topology<Simplex>
{
public:
  TriangulationBase( Points& points , coord_t number ) :
    Topology<Simplex>(points,number,1)
  {}

  virtual void extract() = 0;
};

template<typename type>
class Triangulation : public TriangulationBase
{
public:
  Triangulation( const Topology<type>& topology );

  void extract();

  const Topology<type>& topology() const { return topology_; }

private:
  const Topology<type>& topology_;
  Points points_;
  std::vector<index_t> parents_;
};

} // avro

#endif
