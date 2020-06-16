//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef AVRO_VORONOI_ALGORITHMS_H_
#define AVRO_VORONOI_ALGORITHMS_H_

#include "common/kdtree.h"

#include "element/simplex.h"

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
