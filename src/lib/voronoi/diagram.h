#ifndef AVRO_VORONOI_NEW_DIAGRAM_H_
#define AVRO_VORONOI_NEW_DIAGRAM_H_

#include "common/parallel_for.h"

#include "mesh/field.hpp"

#include "voronoi/cell.h"

#include <geogram/nn_search.h>

namespace avro
{

namespace voronoi
{

class Domain : public Topology<Simplex> {

public:
  Domain( const Topology<Simplex>& topology , coord_t dim ) :
    Topology<Simplex>(lifted_points_,topology.number()),
    lifted_points_(dim)
  {
    // copy the topology
    TopologyBase::copy(topology);

    // copy the points, padding zeros to the additional dimensions if necessary
    avro_assert( dim >= topology.points().dim() );
    std::vector<real_t> x(dim,0.0);
    for (index_t k = 0; k < topology.points().nb(); k++) {
      for (coord_t d = 0; d < topology.points().dim(); d++)
        x[d] = topology.points()[k][d];
      lifted_points_.create(x.data());
    }

    // create the search structures
    search_ = GEO::NearestNeighborSearch::create(dim,"ANN");
    search_->set_points( lifted_points_.nb() , lifted_points_[0] );
    build_structures(); // populates neighbours and inverse
  }

  void get_candidate_elems( const real_t* x , std::vector<index_t>& ball ) {

    // find the closest vertex in the domain to this site
    real_t distance;
    index_t nearest;
    search_->get_nearest_neighbors( 1 , x , &nearest , &distance );
    if (nearest >= lifted_points_.nb()) {
      lifted_points_.print();
      for (coord_t d = 0; d < lifted_points_.dim(); d++)
        printf("x[%u] = %g\n",d,x[d]);
    }
    avro_assert_msg( nearest < lifted_points_.nb() , "nearest = %lu" , nearest );
    inverse().ball( nearest , ball );
  }

private:
  Points lifted_points_;
  GEO::NearestNeighborSearch* search_;
};

class PowerDiagram : public Topology<Polytope> {

public:
  PowerDiagram( const Topology<Simplex>& topology , coord_t dim ) :
    Topology<Polytope>(vertices_,topology.number()),
    sites_(dim),
    vertices_(dim),
    ambient_dim_(topology.number()),
    domain_(topology,dim),
    search_(nullptr),
    verbose_(true)
  {}

  // sets the ambient dimension (2d or 3d)
  void set_ambient_dimension( coord_t dim ) { ambient_dim_ = dim; }
  coord_t ambient_dimension() const { return ambient_dim_; }

  void initialize();

  void allocate_sites( index_t n );
  void set_sites( const real_t* x );
  void set_sites( const Points& p );

  void compute();
  void compute( index_t k );

  bool get_triangles( std::vector<index_t>& triangles , std::vector<index_t>& parent ) const override {
    // overrides get_triangles of Topology since we already know the visualization triangles
    triangles.assign( triangles_.begin() , triangles_.end() );
    parent.assign( triangle2site_.begin() , triangle2site_.end() );
    return true;
  }

  void get_edges( std::vector<index_t>& edges ) const override {
    // overrides get_edges of Topology since we already know the visualization edges
    edges.assign( edges_.begin() , edges_.end() );
  }

  void accumulate(bool complete=false);
  void add_cell( const voronoi::Cell& cell , bool complete );
  void create_field();

  // returns the volume and boundary area so that we can check these in unit tests
  real_t volume() const { return volume_; }
  real_t area() const { return area_; }
  const std::vector<real_t>& cell_volume() const { return cell_volume_; }
  const std::vector<real_t>& centroids() const { return centroid_; }

  // optimizers
  void optimize_points( index_t nb_iter );
  void optimize_points_lloyd( index_t nb_iter );
  void optimize_weights( index_t nb_iter , const std::vector<real_t>& mass );
  bool optimize_weights_kmt( index_t nb_iter , const std::vector<real_t>& mass );

  // optimization-related functions
  void start();
  index_t& iteration() { return iteration_; }
  real_t evaluate_objective( index_t n , const real_t* x , real_t* grad );

  const Points& sites() const { return sites_; }
  Points& sites() { return sites_; }
  index_t nb_cells() const { return cell_.size(); }
  const Cell& cell( index_t k ) const { return *cell_[k].get(); }

  real_t nu( index_t k ) const { return nu_[k]; }
  real_t weight( index_t k ) const { return weight_[k]; }

  real_t regularization() const { return 1e-3; }

  void set_verbose( bool x ) { verbose_ = x; }
  real_t residual() const { return residual_; }

private:

  Points sites_;                // voronoi sites (cell generators)
  Points vertices_;             // voronoi vertices (i.e. where polytopes meet)
  std::vector<real_t> weight_;  // weights imposed on voronoi sites
  coord_t ambient_dim_;         // ambient dimension of the problem (2d or 3d),
                                // note: the sites may be lifted higher if weights are present

  Domain domain_;                           // the domain (a mesh) we will use to clip the voronoi diagram
  std::vector<std::shared_ptr<Cell>> cell_; // voronoi cells (one for each site)
  GEO::NearestNeighborSearch* search_;      // nearest neighbour search structure for the voronoi sites

  real_t volume_;       // volume of the power diagram (used for testing)
  real_t area_;         // area of the boundary of the power diagram (used for testing)

  real_t energy_;                // CVT energy functional
  std::vector<real_t> nu_;       // target mass for the optimal transport problem
  std::vector<real_t> de_dx_;    // deriv. of CVT energy w.r.t. site coordinates
  std::vector<real_t> de_dw_;    // deriv. of CVT energy w.r.t. site weights
  std::vector<real_t> centroid_; // centroids of the power cells
  std::vector<real_t> cell_volume_;

  // optimization-related data
  index_t mode_;         // mode = 0 for coordinate optimiation, mode = 1 for weight optimization
  index_t iteration_;    // objective function evaluation counter
  index_t sub_iteration_;
  real_t time_voronoi_;  // time spent computing the voronoi diagram (seconds)
  bool verbose_;
  real_t residual_;

public:
  // a helper field class to visualize the colour of each voronoi cell
  class SiteField : public Field<Polytope,real_t> {
  public:
    SiteField( PowerDiagram& diagram ) :
      Field<Polytope,real_t>(diagram,0,DISCONTINUOUS)
    {
      build();
      for (index_t k = 0; k < diagram.nb(); k++)
        this->value(k) = real_t(diagram.polytope2site_[k]);
    }

    index_t nb_rank() const { return 1; }
  };

  SiteField& site_field() { return *site_field_.get(); }
  const std::vector<index_t>& polytope2site() const { return polytope2site_; }

private:
  // visualization-related data
  std::vector<index_t> triangles_;
  std::vector<index_t> triangle2site_;
  std::vector<index_t> edges_;
  std::vector<index_t> polytope2site_;
  std::shared_ptr<SiteField> site_field_;
};

} // voronoi

} // avro

#endif
