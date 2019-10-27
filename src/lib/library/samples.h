#ifndef URSA_LIB_LIBRARY_SAMPLES_H_
#define URSA_LIB_LIBRARY_SAMPLES_H_

#include "mesh/topology.h"
#include "mesh/points.h"

namespace ursa
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

} // ursa

#endif
