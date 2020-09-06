#ifndef AVRO_LIB_PHYSICS_PARTICLES_H_
#define AVRO_LIB_PHYSICS_PARTICLES_H_

#include "voronoi/power.h"

namespace avro
{

class PowerFacets;

class ParticleFacets : public Topology<Polytope>
{
public:
  ParticleFacets( Points& vertices , const PowerFacets& facets );

  void extract( const PowerFacets& facets );

private:
  std::vector<int> cellL_;
  std::vector<int> cellR_;

};

class Particles : public Topology<Polytope>
{
public:
  Particles( Points& points );

  void extract();

private:
  Points& points_;
  Points vertices_;
  ParticleFacets facets_;

  delaunay::PowerDiagram power_diagram_;
};

#endif
