#ifndef avro_LIB_LIBRARY_SAMPLES_H_
#define avro_LIB_LIBRARY_SAMPLES_H_

#include "mesh/topology.h"
#include "mesh/points.h"

namespace avro
{

namespace library
{

class TwoTriangles : public Topology<Simplex>
{
public:
  TwoTriangles();

  Topology<Simplex>& edges() { return edges_; }

private:
  Points points_;
  Topology<Simplex> edges_;
};

} // library

} // avro

#endif
