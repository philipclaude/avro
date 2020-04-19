#ifndef AVRO_VORONOI_ALGORITHMS_H_
#define AVRO_VORONOI_ALGORITHMS_H_

#include "master/simplex.h"

#include "mesh/topology.h"

namespace avro
{


class BowerWatson : public Topology<Simplex>
{
public:
  BowerWatson( Points& points );

  void compute();

  void initialize();

  int insphere( index_t elem , index_t point );

private:
  Points points_;
  Points& delaunay_;

};

}

#endif
