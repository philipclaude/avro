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
    if (nearest >= nb()) {
      lifted_points_.print();
      for (coord_t d = 0; d < lifted_points_.dim(); d++)
        printf("x[%u] = %g\n",d,x[d]);
    }
    avro_assert( nearest < nb() );
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
    search_(nullptr)
  {}

  // sets the ambient dimension (2d or 3d)
  void set_ambient_dimension( coord_t dim ) { ambient_dim_ = dim; }

  // creates voronoi cells and sites the target mass and weight arrays
  // should only be called once
  void initialize() {
    cell_.resize( sites_.nb() );
    for (index_t k = 0; k < sites_.nb(); k++) {
      cell_[k] = std::make_shared<Cell>( k , sites_ , domain_ , *search_ );
      cell_[k]->set_ambient_dimension(ambient_dim_);
    }
    nu_.resize(sites_.nb() , 0.0 );
    weight_.resize( sites_.nb() , 0.0 );
  }

  void set_sites( const real_t* x ) {

    // copy the points
    index_t i = 0;
    for (index_t k = 0; k < sites_.nb(); k++)
    for (coord_t d = 0; d < ambient_dim_; d++)
      sites_[k][d] = x[i++];

    // set the points into the nearest neighbour search structure
    if (search_ == nullptr)
      search_ = GEO::NearestNeighborSearch::create(sites_.dim(),"ANN");
    search_->set_points( sites_.nb() , sites_[0] );
  }

  void set_sites( const Points& p ) {

    coord_t dim = p.dim();
    avro_assert( dim == sites_.dim() );

    // copy the points and set the coordinates into the search structure
    p.copy(sites_);
    if (search_ == nullptr)
      search_ = GEO::NearestNeighborSearch::create(dim,"ANN");
    search_->set_points( sites_.nb() , sites_[0] );
  }

  // set the site weights (lifts the sites to ambient_dim+1)
  void set_weights( const real_t* x );

  void compute() {

    // clear any previous power diagram data
    Topology<Polytope>::clear();
    vertices_.clear();
    triangles_.clear();
    triangle2site_.clear();
    edges_.clear();
    polytope2site_.clear();

    // clear any optimization data
    de_dx_.resize( sites_.nb() * ambient_dim_ , 0.0 );
    de_dw_.resize( sites_.nb() , 0.0 );
    centroid_.resize( sites_.nb() * ambient_dim_ , 0.0 );

    // reset the energy, volume and area
    energy_ = 0.0;
    volume_ = 0.0;
    area_   = 0.0;

    // compute the power diagra
    clock_t t0 = clock();
    #if 0
    typedef PowerDiagram thisclass;
    ProcessCPU::parallel_for(
      parallel_for_member_callback( this , &thisclass::compute ), 0,cell_.size()
    );
    #else
    for (index_t k = 0; k < cell_.size(); k++) {
      compute(k);
    }
    #endif

    // record the time it took to compute the voronoi diagram
    time_voronoi_ = real_t(clock()-t0)/real_t(CLOCKS_PER_SEC)/ProcessCPU::maximum_concurrent_threads();
  }

  void compute( index_t k ) {

    // get the candidate list of elements in the background mesh to clip against
    // (only one will be used)
    std::vector<index_t> ball;
    domain_.get_candidate_elems( sites_[k] , ball );

    // clip the voronoi cell against the mesh
    cell_[k]->compute( ball );

    // retrieve the mass and moment data
    const std::vector<real_t>& moment = cell_[k]->moment();
    real_t volume_k = cell_[k]->volume();

    #if 0 // this check should only be used when weights are 0, oterwise cells could indeed vanish
    if (volume_k <= 0.0) {
      sites_.print(k);
    }
    avro_assert( volume_k > 0.0 );
    #endif

    // accumulate the total volume and boundary area (used for testing)
    volume_ += volume_k;
    area_   += cell_[k]->boundary_area();

    // compute the centroid of the cell
    if (volume_k > 0.0) {
      for (coord_t d = 0; d < ambient_dim_; d++)
        centroid_[k*ambient_dim_+d] = moment[d]/volume_k;
    }
    else {
      for (coord_t d = 0; d < ambient_dim_; d++)
        centroid_[k*ambient_dim_+d] = sites_[k][d]; // retain the previous value
    }

    // compute the gradient of the energy with respect to the site locations
    for (coord_t d = 0; d < ambient_dim_; d++)
      de_dx_[k*ambient_dim_+d] = 2.0*( volume_k*sites_[k][d] - moment[d] );

    // compute the gradient of the energy with respect to the weights
    de_dw_[k] = volume_k - nu_[k];

    // add the contribution of the cell energy to the total
    energy_   += cell_[k]->energy() + weight_[k] * (volume_k - nu_[k]);
  }

  void accumulate() {
    // accumulate all the cells into a polytope mesh that we can visualize
    for (index_t k = 0; k < cell_.size(); k++) {
      add_cell( *cell_[k].get() );
      cell_[k]->clear();
    }
  }

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

  void add_cell( const voronoi::Cell& cell ) {

    // add the points
    index_t nb_points = vertices_.nb();
    for (index_t j = 0; j < cell.points().nb(); j++) {
      vertices_.create( cell.points()[j] );
    }

    // add a dummy polytope
    // (we don't actually need the polytope since we already have the triangles to visualize)
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

  void create_field() {
    // creates a field so that we can visualize each voronoi cell with a constant colour
    site_field_ = std::make_shared<SiteField>(*this);
    fields().make("sites",site_field_);
  }

  // returns the volume and boundary area so that we can check these in unit tests
  real_t volume() const { return volume_; }
  real_t area() const { return area_; }

  // optimizers
  void optimize_points( index_t nb_iter );
  void optimize_points_lloyd( index_t nb_iter );
  void optimize_weights( index_t nb_iter , const std::vector<real_t>& mass );

  // optimization-related functions
  void start();
  index_t& iteration() { return iteration_; }
  real_t evaluate_objective( index_t n , const real_t* x , real_t* grad );

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

  // optimization-related data
  index_t mode_;         // mode = 0 for coordinate optimiation, mode = 1 for weight optimization
  index_t iteration_;    // objective function evaluation counter
  real_t time_voronoi_;  // time spent computing the voronoi diagram (seconds)

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
