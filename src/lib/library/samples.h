#ifndef luma_LIB_LIBRARY_SAMPLES_H_
#define luma_LIB_LIBRARY_SAMPLES_H_

#include "mesh/topology.h"
#include "mesh/points.h"

namespace luma
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

} // luma

#endif
