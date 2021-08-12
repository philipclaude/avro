#ifndef AVRO_VORONOI_NEW_DIAGRAM_H_
#define AVRO_VORONOI_NEW_DIAGRAM_H_

#include "common/parallel_for.h"

#include "mesh/field.hpp"

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
    search_(nullptr),
    volume_(0.0),
    boundary_area_(0.0),
    ambient_dim_(topology.number())
  {}

  void set_ambient_dimension( coord_t dim ) { ambient_dim_ = dim; }

  void initialize() {
    cell_.resize( sites_.nb() );
    for (index_t k = 0; k < sites_.nb(); k++) {
      cell_[k] = std::make_shared<Cell>( k , sites_ , domain_ , *search_ );
      cell_[k]->set_ambient_dimension(ambient_dim_);
    }
    nu_.resize(sites_.nb() , 0.0 );
    weights_.resize( sites_.nb() , 0.0 );
  }

  void set_sites( const real_t* x ) {
    index_t i = 0;
    for (index_t k = 0; k < sites_.nb(); k++)
    for (coord_t d = 0; d < sites_.dim(); d++)
      sites_[k][d] = x[i++];

    if (search_ == nullptr)
      search_ = GEO::NearestNeighborSearch::create(sites_.dim(),"ANN");
    search_->set_points( sites_.nb() , sites_[0] );
  }

  void set_sites( const Points& p ) {

    coord_t dim = p.dim();
    avro_assert( dim == sites_.dim() );

    // copy the points and set the search structure
    p.copy(sites_);
    if (search_ == nullptr)
      search_ = GEO::NearestNeighborSearch::create(dim,"ANN");
    search_->set_points( sites_.nb() , sites_[0] );
  }
  void set_weights( const real_t* x );

  void set_target_mass( const std::vector<real_t>& nu ) { nu_ = nu; }

  void compute() {

    Topology<Polytope>::clear();
    vertices_.clear();
    triangles_.clear();
    triangle2site_.clear();
    edges_.clear();
    polytope2site_.clear();

    clock_t t0 = clock();

    de_dx_.resize( sites_.nb() * ambient_dim_ , 0.0 );
    de_dw_.resize( sites_.nb() , 0.0 );
    energy_ = 0.0;

    #if 1
    typedef PowerDiagram thisclass;
    ProcessCPU::parallel_for(
      parallel_for_member_callback( this , &thisclass::compute ), 0,cell_.size()
    );
    #else
    for (index_t k = 0; k < cell_.size(); k++) {
      compute(k);
    }
    #endif

    clock_t t1 = clock();
    time_voronoi_ = real_t(t1-t0)/real_t(CLOCKS_PER_SEC)/ProcessCPU::maximum_concurrent_threads();
  }

  void compute( index_t k ) {
    std::vector<index_t> ball;
    domain_.get_candidate_elems( sites_[k] , ball );
    cell_[k]->compute( ball );

    const std::vector<real_t>& moment = cell_[k]->moment();
    real_t mass = cell_[k]->volume();

    for (coord_t d = 0; d < ambient_dim_; d++)
      de_dx_[k*ambient_dim_+d] = 2.0*( mass*sites_[k][d] - moment[d] );

    de_dw_[k] = mass - nu_[k];
    energy_ += cell_[k]->energy() + weights_[k] * (mass - nu_[k]);
  }

  void accumulate() {
    for (index_t k = 0; k < cell_.size(); k++) {
      add_cell( *cell_[k].get() );
      cell_[k]->clear();
    }
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
    boundary_area_ += cell.boundary_area();

    index_t nb_points = vertices_.nb();
    for (index_t j = 0; j < cell.points().nb(); j++) {
      vertices_.create( cell.points()[j] );
    }

    // add the polytope
    for (index_t j = 0; j < cell.nb(); j++) {
      index_t dummy = 0;
      add( &dummy , 1 );
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
  real_t boundary_area() const { return boundary_area_; }

  void start();
  void optimize_points( index_t nb_iter );
  void optimize_weights( index_t nb_iter , const std::vector<real_t>& mass );

  index_t& iteration() { return iteration_; }

  real_t evaluate_objective( index_t n , const real_t* x , real_t* grad );

private:

  Points sites_;
  Points vertices_;
  std::vector<index_t> triangles_; // for visualization
  std::vector<index_t> triangle2site_;
  std::vector<index_t> edges_;

  Domain domain_;
  std::vector<std::shared_ptr<Cell>> cell_;
  GEO::NearestNeighborSearch* search_;

  real_t volume_;
  real_t boundary_area_;
  coord_t ambient_dim_;

  real_t energy_;
  std::vector<real_t> nu_; // target mass
  std::vector<real_t> weights_;
  std::vector<real_t> de_dx_;
  std::vector<real_t> de_dw_;
  index_t mode_;
  index_t iteration_;
  real_t time_voronoi_;

  class SiteField : public Field<Polytope,real_t> {
  public:
    SiteField( PowerDiagram& diagram ) :
      Field<Polytope,real_t>(diagram,0,DISCONTINUOUS)
    {
      build();
      for (index_t k = 0; k < diagram.nb(); k++)
        this->value(k) = real_t(diagram.polytope2site(k));
      //this->dof().print();
    }

    index_t nb_rank() const { return 1; }
  };

  std::vector<index_t> polytope2site_;
  std::shared_ptr<SiteField> site_field_;
};

} // voronoi

} // avro

#endif
