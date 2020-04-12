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
  void generate_sites() { avro_implement; }

private:
  Points points_;  // to accumulate the points in the calculation
  Delaunay sites_;

  bool hierarchical_;
  bool exact_;

  std::vector<std::shared_ptr<Topology<Simplex>>> topologies_;
  std::vector<Entity*> entities_;
  std::vector<std::shared_ptr<VoronoiSites>> sites_fields_;
};

} // delaunay

} // avro

#endif
