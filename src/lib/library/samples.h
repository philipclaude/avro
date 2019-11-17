#ifndef LUNA_LIB_LIBRARY_SAMPLES_H_
#define LUNA_LIB_LIBRARY_SAMPLES_H_

#include "mesh/topology.h"
#include "mesh/points.h"

namespace luna
{

namespace library
{

class TwoTriangles : public Topology<Simplex<Lagrange>>
{
public:
  TwoTriangles();

  Topology<Simplex<Lagrange>>& edges() { return edges_; }

private:
  Points points_;
  Topology<Simplex<Lagrange>> edges_;
};

} // library

} // luna

#endif
