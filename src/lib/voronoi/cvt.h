//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef AVRO_LIB_VORONOI_CVT_H_
#define AVRO_LIB_VORONOI_CVT_H_

#include "mesh/points.h"
#include "mesh/topology.h"

#include "voronoi/delaunay.h"

namespace avro
{

namespace delaunay
{

class RestrictedVoronoiDiagram;
class VoronoiSites;

class CentroidalVoronoiTessellation : public Topology<Polytope>
{
public:
  CentroidalVoronoiTessellation( const Topology<Simplex>& topology , Points& sites , bool hierarchical=false );

  void compute( index_t nb_iter );
  void generate_sites( int nb_samples=-1 );

private:
  Points points_;  // to accumulate the points in the calculation
  Delaunay sites_;

  bool hierarchical_;
  bool exact_;

  std::vector<std::shared_ptr<Topology<Simplex>>> topologies_;
  std::vector<Entity*> entities_;
  std::vector<std::shared_ptr<VoronoiSites>> sites_fields_;

  void sample_geometry( Entity* entity , index_t nb_samples );
};

} // delaunay

} // avro

#endif
