#ifndef AVRO_LIB_VORONOI_CELL_H_
#define AVRO_LIB_VORONOI_CELL_H_

#include "element/polytope.h"

#include "mesh/topology.h"

#include "voronoi/voronoi.h"

#include <memory>

#include <geogram/nn_search.h>

namespace avro
{

namespace delaunay
{

class RVDFacets;

class VoronoiCell : public Topology<Polytope>
{
public:
  VoronoiCell( index_t site , const Delaunay& delaunay , const NearestNeighbours& neighbours ,
               const TopologyBase& domin , bool exact , bool simplex ,
               GEO::NearestNeighborSearch* nns=nullptr , index_t nb_nns=0 );

  void initialize( const std::vector<index_t>& E = std::vector<index_t>() );
  void compute(const std::vector<index_t>& E = std::vector<index_t>());

  void set_facets( RVDFacets* facets ) { facets_ = facets; }

  // for passing initial guess back to neighbour reconstruction
  index_t nb_neighbours() const { return nn_.size(); }
  const std::vector<index_t>& neighbours() const { return nn_; }

  void generate_simplices();
  const Topology<Simplex>& simplices() const { return simplices_; }

  real_t time_neighbours() const { return time_neighbours_; }
  real_t time_clip() const { return time_clip_; }
  real_t time_decompose() const { return time_decompose_; }

  bool incomplete() const { return incomplete_; }

  void get_bisector( int b , index_t& p0 , index_t& p1 ) const;

  void print() const;

private:
  void initialize_polytope( const std::vector<index_t>& E );
  void initialize_simplex();
  void compute_polytope();
  void compute_simplex();

  void enlarge_neighbours();

  void clip_by_bisector( index_t j , index_t bj );
  int  clip_edge( index_t e0 , index_t e1 , const int b , std::vector<index_t>& q ,int& q0, int& q1  );
  bool security_radius_reached( index_t bj ) const;

  int add_bisector( index_t p0 , index_t p1 );

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

  GEO::NearestNeighborSearch* nns_;
  std::vector<index_t> nn_;
  bool recycle_neighbours_;

  Points simplex_points_;
  Topology<Simplex> simplices_;

  real_t time_neighbours_;
  real_t time_clip_;
  real_t time_decompose_;

  std::vector<index_t> qpolytope_;

  bool incomplete_;

  std::vector<index_t> pedges_;
  std::vector<index_t> qedges_;
  std::vector<index_t> qplane_;
};

class VoronoiDiagram : public Topology<Polytope>
{
  typedef VoronoiDiagram thisclass;

public:
  VoronoiDiagram( Delaunay& delaunay , const TopologyBase& domain , bool simplex=false );

  void compute( bool exact );
  void clip( const index_t k )
  {
    cells_[k]->compute(domain_edges_);
  }

  const std::vector<index_t>& sites() const { return sites_; }
  const VoronoiCell& cell( index_t k ) const { return *cells_[k].get(); }

  index_t vertex2site( index_t k ) const { return vertex2site_[k]; }
  const std::vector<SymbolicVertex>& symbolic_vertices() const { return symbolic_vertices_; }
  const std::map<int,Bisector>& bisectors() const { return bisectors_; }

  real_t time_neighbours() const { return time_neighbours_; }
  real_t time_voronoi() const { return time_voronoi_; }

private:
  Points points_;
  Delaunay& delaunay_;
  const TopologyBase& domain_;
  std::vector<std::shared_ptr<VoronoiCell>> cells_;
  std::vector<index_t> sites_;
  bool simplex_;
  real_t time_neighbours_;
  real_t time_voronoi_;

  GEO::NearestNeighborSearch* nns_;

  std::vector<index_t> vertex2site_;
  std::vector<SymbolicVertex> symbolic_vertices_;
  std::map<int,Bisector> bisectors_;

  std::vector<index_t> domain_edges_;

};

} // delaunay

} // avro

#endif
