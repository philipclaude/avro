#ifndef AVRO_LIB_VORONOI_CELL_H_
#define AVRO_LIB_VORONOI_CELL_H_

#include "element/polytope.h"

#include "mesh/topology.h"

#include "voronoi/voronoi.h"

#include <memory>

namespace avro
{

namespace delaunay
{

class RVDFacets;

class VoronoiCell : public Topology<Polytope>
{
public:
  VoronoiCell( index_t site , const Delaunay& delaunay , const NearestNeighbours& neighbours , const TopologyBase& domin , bool exact , bool simplex );

  void initialize();
  void compute();

  void set_facets( RVDFacets* facets ) { facets_ = facets; }

private:
  void initialize_polytope();
  void initialize_simplex();
  void compute_polytope();
  void compute_simplex();

  void clip_by_bisector( index_t j , index_t bj );
  void clip_edge( index_t e0 , index_t e1 , const int b , std::vector<index_t>& q );
  bool security_radius_reached( index_t bj ) const;

  int add_bisector( index_t p0 , index_t p1 );
  void get_bisector( int b , index_t& p0 , index_t& p1 ) const;

  const index_t site_;
  Points points_;
  const Delaunay& delaunay_;
  const NearestNeighbours& neighbours_;
  const TopologyBase& domain_;
  const bool exact_;

  std::vector<Vertex> vertex_;
  std::vector<index_t> polytope_;  // current polytope points

  RVDFacets* facets_;
  bool simplex_;

  std::map<Bisector,int> bisector_;
  std::map<int,Bisector> ids_;
};

class VoronoiDiagram : public Topology<Polytope>
{
  typedef VoronoiDiagram thisclass;

public:
  VoronoiDiagram( Delaunay& delaunay , const TopologyBase& domain , bool simplex=false );

  void compute( bool exact );
  void clip( const index_t k )
    { cells_[k]->compute(); }

  const std::vector<index_t>& sites() const { return sites_; }

private:
  Points points_;
  Delaunay& delaunay_;
  const TopologyBase& domain_;
  std::vector<std::shared_ptr<VoronoiCell>> cells_;
  std::vector<index_t> sites_;
  bool simplex_;
};

} // delaunay

} // avro

#endif
