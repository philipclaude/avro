#ifndef avro_LIB_VORONOI_GEOMETRY_H_
#define avro_LIB_VORONOI_GEOMETRY_H_

#include "element/polytope.h"

#include "mesh/points.h"
#include "mesh/topology.h"

#include <memory>
#include <vector>

namespace avro
{

class Entity;
class Delaunay;
template<typename type> class Boundary;

namespace delaunay
{
class RestrictedVoronoiDiagram;
}

class GeometryConformingRVD : public Topology<Polytope>
{
public:

  typedef Topology<Simplex> RestrictedDelaunayTriangulation;

  GeometryConformingRVD( const Topology<Simplex>& topology , Delaunay& sites );

  void initialize();
  void compute();
  void extract_triangulations();
  Topology<Simplex>& triangulation() { return triangulation_; }

private:
  const Topology<Simplex>& topology_;
  Delaunay& sites_;
  Points vertices_;
  Topology<Simplex> triangulation_;
  std::shared_ptr< Boundary<Simplex> > boundary_;

  std::vector< std::shared_ptr< delaunay::RestrictedVoronoiDiagram > > rvd_;
  std::vector< std::shared_ptr< RestrictedDelaunayTriangulation > > rdt_;
  std::vector< Entity* > entity_;
};

} // avro

#endif
