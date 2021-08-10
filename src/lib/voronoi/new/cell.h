#ifndef AVRO_LIB_VORONOI_POWER_CELL_H_
#define AVRO_LIB_VORONOI_POWER_CELL_H_

#include "element/polytope.h"
#include "mesh/topology.h"

#include "voronoi/new/facets.h"
#include "voronoi/new/vertex.h"

#include <nnsearch/nn_search.h>

namespace avro
{

template<typename type> class Topology;
class Simplex;
class Parameters;

namespace voronoi
{

class Cell : public Topology<Polytope> {

public:
  Cell( index_t site , const Points& sites ,
        const Topology<Simplex>& domain , GEO::NearestNeighborSearch& searcher );

  void set_parameters( const Parameters& params );
  void compute( index_t elem );

  real_t get_mass() const;       // for gradient
  real_t get_moment() const;     // for gradient
  real_t get_face_terms() const; // for hessian

  // for plotting the resulting power diagram
  //void get_points( std::vector<real_t>& points ) const;
  //void get_edges( std::vector<index_t>& edges ) const;
  //void get_triangles( std::vector<index_t>& triangles ) const;

protected:

  void initialize( index_t elem );
  void enlarge_neighbours();

  void clip( index_t elem );
  void clip_simplex( index_t elem );
  void clip_by_bisector( index_t j , index_t bj );
  int  clip_edge( index_t e0 , index_t e1 , const int b , std::vector<index_t>& q ,int& q0, int& q1  );
  bool security_radius_reached( index_t bj ) const;

  void get_bisector( int b , index_t& p0 , index_t& p1 ) const;
  int  add_bisector( index_t p0 , index_t p1 );

  void generate_simplices();

  template<typename Integrand>
  void integrate( const Integrand& integrand );

private:
  index_t site_; // the index of the voronoi site we are computing the cell of (gl_PrimitiveID)
  const Points& delaunay_;             // the delaunay vertices/voronoi sites (texture)
  const Topology<Simplex>& domain_; // topology from which we will grab an element to clip against (texture)
  GEO::NearestNeighborSearch& search_; // search structure through delaunay points that may need to be enlarged (texture)
  //const Facets& facets_; // mesh facets (in global numbering) (texture)

  std::vector<index_t> neighbours_;     // current list of nearest neighbours to the site
  std::unordered_set<index_t> visited_; // set of elements which have been visited
  std::map<Bisector,int> bisector_;     // map from actual bisector to bisector label
  std::map<int,Bisector> ids_;          // map from bisector label to actual bisector

  std::vector<Vertex>  vertex_;    // list of vertices
  std::vector<index_t> polytope_;  // current polytope points
  std::vector<index_t> qpolytope_; // temporary polytope used while clipping
  std::vector<index_t> pedges_;    // current polytope edges
  std::vector<index_t> qedges_;    // temporary polytopes edges used while clipping
  std::vector<index_t> qplane_;    // points on the hyperplane of the current clipped polytope

  Points points_;
  bool exact_;
};

} // voronoi

} // avro

#endif
