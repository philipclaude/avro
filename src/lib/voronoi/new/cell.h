#ifndef AVRO_LIB_VORONOI_POWER_CELL_H_
#define AVRO_LIB_VORONOI_POWER_CELL_H_

#include "element/polytope.h"
#include "mesh/topology.h"

#include "voronoi/new/facets.h"
#include "voronoi/new/vertex.h"

#include <geogram/nn_search.h>

namespace avro
{

template<typename type> class Topology;
class Simplex;
class Parameters;

namespace voronoi
{

typedef struct {

  index_t elem; // background mesh element from which the facet originated

  int pi; // nonnegative
  int pj; // could be negative if a boundary face

  std::vector<real_t> qi; // site coordinates
  std::vector<real_t> bi; // cell centroid
  std::vector<real_t> Abij; // cell moment (Aij * bij)

  real_t Aij; // area of the face
  real_t dij; // distance from qi to face
  real_t lij; // distance from qi to qj
  std::vector<real_t> bij; // centroid of the face shared with site j
  std::vector<real_t> normal; // normal to the facet

} PowerFacet;

class Cell : public Topology<Polytope> {

public:
  Cell( index_t site , const Points& sites ,
        const Topology<Simplex>& domain , GEO::NearestNeighborSearch& searcher );

  void set_parameters( const Parameters& params );
  void compute( const std::vector<index_t>& candidates );
  void set_ambient_dimension( coord_t dim ) { ambient_dim_ = dim; }

  real_t get_mass() const;       // for gradient
  real_t get_moment() const;     // for gradient
  real_t get_face_terms() const; // for hessian

  const std::vector<index_t>& triangles() const { return triangles_; }
  const std::vector<index_t>& edges() const { return edges_; }

  index_t site() const { return site_; }
  real_t volume() const { return volume_; }

  void clear();

  index_t decode_mesh_facet( index_t elem , int b ) {
    index_t id = -b - 1;
    avro_assert( id / (domain_.number()+1) == elem );
    index_t f = id % (domain_.number()+1);
    return f;
  }

  int encode_mesh_facet( index_t f ) {
    index_t id = elem_ * (domain_.number()+1) + f;
    return -id -1;
  }

  void finish_facets();
  real_t boundary_area() const { return boundary_area_; }
  const std::vector<real_t>& moment() { return moment_; }
  real_t compute_energy( const std::vector<const real_t*>& X ) const;
  real_t energy() const { return energy_; }

  const std::map<int,PowerFacet>& facets() const { return facets_; }

protected:

  void initialize( index_t elem );
  void enlarge_neighbours();

  bool clip_simplex( index_t elem );
  void clip_by_bisector( index_t j , index_t bj );
  int  clip_edge( index_t e0 , index_t e1 , const int b , std::vector<index_t>& q ,int& q0, int& q1  );
  bool security_radius_reached( index_t bj ) const;

  void get_bisector( int b , index_t& p0 , index_t& p1 ) const;
  int  add_bisector( index_t p0 , index_t p1 );

  void generate_simplices();

private:
  index_t site_; // the index of the voronoi site we are computing the cell of (gl_PrimitiveID)
  const Points& delaunay_;             // the delaunay vertices/voronoi sites (texture)
  const Topology<Simplex>& domain_; // topology from which we will grab an element to clip against (texture)
  GEO::NearestNeighborSearch& search_; // search structure through delaunay points that may need to be enlarged (texture)

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

  // visualization data
  std::vector<index_t> triangles_;
  std::vector<index_t> edges_;

  index_t elem_;  // current mesh elements against which we are clipping
  Points points_; // points used to describe the polytopes
  bool exact_;    // whether to use exact arithmetic
  coord_t ambient_dim_; // ambient dimension of the mesh (either 2d or 3d)

  // cell terms
  real_t volume_;
  std::vector<real_t> moment_;
  real_t energy_;

  // face terms
  std::map<int,PowerFacet> facets_;
  real_t boundary_area_; // this is only used for unit testing
};

} // voronoi

} // avro

#endif
