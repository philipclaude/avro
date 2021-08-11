#ifndef AVRO_VORONOI_NEW_DIAGRAM_H_
#define AVRO_VORONOI_NEW_DIAGRAM_H_

#include "voronoi/new/cell.h"

#include <nnsearch/nn_search.h>

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
    std::vector<real_t> X( topology.points().nb()*dim , 0.0 );
    for (index_t k = 0; k < topology.points().nb(); k++) {
      for (coord_t d = 0; d < topology.points().dim(); d++) {
        x[d] = topology.points()[k][d];
        X[k*dim+d] = x[d];
      }
      lifted_points_.create(x.data());
    }

    // create the search structures
    search_ = GEO::NearestNeighborSearch::create(dim,"ANN");
    search_->set_points( lifted_points_.nb() , X.data() );
    build_structures(); // populates neighbours and inverse
  }

  void get_candidate_elems( const real_t* x , std::vector<index_t>& ball ) {

    // find the closest vertex in the domain to this site
    real_t distance;
    index_t nearest;
    search_->get_nearest_neighbors( 1 , x , &nearest , &distance );
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
    domain_(topology,dim),
    volume_(0.0)
  {
    search_ = GEO::NearestNeighborSearch::create(dim,"ANN");
  }

  void initialize( index_t nb_cells ) {
    cell_.resize( nb_cells );
    for (index_t k = 0; k < nb_cells; k++)
      cell_[k] = std::make_shared<Cell>( k , sites_ , domain_ , *search_ );
  }

  void set_sites( const real_t* x );
  void set_sites( const Points& p ) {

    p.copy(sites_);

    coord_t dim = p.dim();
    avro_assert( dim == sites_.dim() );
    std::vector<real_t> X( sites_.nb()*dim , 0.0 );
    for (index_t k = 0; k < sites_.nb(); k++) {
      for (coord_t d = 0; d < dim; d++) {
        X[k*dim+d] = sites_[k][d];
      }
    }
    search_->set_points( sites_.nb() , X.data() );
  }
  void set_weights( const real_t* x );

  void compute() {
    for (index_t k = 0; k < cell_.size(); k++) {
      compute(k);
    }
  }

  void compute( index_t k ) {
    std::vector<index_t> ball;
    domain_.get_candidate_elems( sites_[k] , ball );
    cell_[k]->compute( ball );
  }

  void accumulate() {
    for (index_t k = 0; k < cell_.size(); k++)
      add_cell( *cell_[k].get() );
  }

  bool get_triangles( std::vector<index_t>& triangles , std::vector<index_t>& parent ) const override {
    triangles.assign( triangles_.begin() , triangles_.end() );
    parent.assign( triangle2site_.begin() , triangle2site_.end() );
    return true;
  }
  void get_edges( std::vector<index_t>& edges ) const override {
    edges.assign( edges_.begin() , edges_.end() );
  }

  void add_cell( const voronoi::Cell& cell ) {

    volume_ += cell.volume();

    index_t nb_points = vertices_.nb();
    for (index_t j = 0; j < cell.points().nb(); j++) {
      vertices_.create( cell.points()[j] );
    }

    // add the polytope
    for (index_t j = 0; j < cell.nb(); j++) {
      std::vector<index_t> polytope = cell.get(j);
      for (index_t i = 0; i < polytope.size(); i++)
        polytope[i] += nb_points;
      add( polytope.data() , polytope.size() );
      polytope2site_.push_back(cell.site());
    }

    // add the triangle data
    const std::vector<index_t>& t = cell.triangles();
    for (index_t j = 0; j < t.size(); j++)
      triangles_.push_back(t[j]+nb_points);

    for (index_t j = 0; j < t.size()/3; j++)
      triangle2site_.push_back(cell.site());

    // add the edge data
    const std::vector<index_t>& e = cell.edges();
    for (index_t j = 0; j < e.size(); j++)
      edges_.push_back(e[j]+nb_points);
  }

  index_t polytope2site( index_t k ) const { return polytope2site_[k]; }

  void create_field() {
    site_field_ = std::make_shared<SiteField>(*this);
    fields().make("sites",site_field_);
  }

  real_t volume() const { return volume_; }

private:

  Points sites_;
  Points vertices_;
  std::vector<index_t> triangles_; // for visualization
  std::vector<index_t> triangle2site_;
  std::vector<index_t> edges_;

  real_t volume_;

  Domain domain_;
  std::vector<std::shared_ptr<Cell>> cell_;
  GEO::NearestNeighborSearch* search_;

  class SiteField : public Field<Polytope,real_t> {
  public:
    SiteField( PowerDiagram& diagram ) :
      Field<Polytope,real_t>(diagram,0,DISCONTINUOUS)
    {
      build();
      for (index_t k = 0; k < diagram.nb(); k++)
        this->value(k) = real_t(diagram.polytope2site(k));
    }

    index_t nb_rank() const { return 1; }
  };

  std::vector<index_t> polytope2site_;
  std::shared_ptr<SiteField> site_field_;
};

} // voronoi

} // avro

#endif
