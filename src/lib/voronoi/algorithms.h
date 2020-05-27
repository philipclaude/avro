#ifndef AVRO_VORONOI_ALGORITHMS_H_
#define AVRO_VORONOI_ALGORITHMS_H_

#include "common/kdtree.h"

#include "master/simplex.h"

#include "mesh/topology.h"

namespace avro
{

class BowyerWatson : public Topology<Simplex>
{
public:
  BowyerWatson( Points& points );

  void compute();

  void initialize();

  real_t insphere( index_t elem , index_t point );

  void search( index_t point , index_t elem , std::set<index_t>& elems );

private:
  Points points_;
  Points& delaunay_;

  std::vector<index_t> fake_;

  PointCloud cloud_;
  std::shared_ptr<KdTreeNd> kdtree_;
};

}

#endif
