#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"

#include "mesh/field.hpp"

#include "voronoi/new/cell.h"

#include <nnsearch/nn_search.h>

using namespace avro;

UT_TEST_SUITE( voronoi_cell_test_suite )

class PowerDiagram : public Topology<Polytope> {

public:
  PowerDiagram( coord_t number , coord_t dim ) :
    Topology<Polytope>(vertices_,number),
    vertices_(dim),
    volume_(0.0)
  {}

  ~PowerDiagram() {}

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
  Points vertices_;
  std::vector<index_t> triangles_; // for visualization
  std::vector<index_t> triangle2site_;
  std::vector<index_t> edges_;

  real_t volume_;

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

UT_TEST_CASE( test_2d )
{
  GEO::PCK::initialize();

  static coord_t number = 2;
  static coord_t dim = number;
  index_t nb_points = 1e3;

  std::vector<index_t> dims(number,10);
  CKF_Triangulation ckf(dims);

  Points points(dim);
  real_t x0[] = {0,0,0}; points.create(x0);
  real_t x1[] = {1,0,0}; points.create(x1);
  real_t x2[] = {0,1,0}; points.create(x2);
  if (number > 2) {
    real_t x3[] = {0,0,1}; points.create(x3);
  }

  #if 0
  Topology<Simplex> topology(points,number);
  index_t t[] = {0,1,2,3};
  topology.add(t,number+1);
  #else
  Topology<Simplex>& topology = ckf;
  #endif
  topology.build_structures();

  Points sites(dim);

  #if 0
  ckf.points().copy( sites );
  nb_points = sites.nb();
  #else
  for (index_t k = 0; k < nb_points; k++) {
    for (coord_t d = 0; d < dim; d++)
      x0[d] = random_within(0.0,1.0);
    sites.create(x0);
  }
  #endif

  GEO::NearestNeighborSearch* searcher = GEO::NearestNeighborSearch::create(dim,"ANN");
  std::vector<real_t> X(sites.nb()*dim);
  for (index_t k = 0; k < sites.nb(); k++) {
    for (index_t d = 0; d < dim; d++)
      X[k*dim+d] = sites[k][d];
  }
  searcher->set_points( sites.nb() , X.data() );

  PowerDiagram diagram(number,dim);

  // create a search tree for the domain points
  GEO::NearestNeighborSearch* domain_searcher = GEO::NearestNeighborSearch::create(dim,"ANN");
  std::vector<real_t> Y( topology.points().nb()*dim );
  for (index_t k = 0; k < topology.points().nb(); k++)
  for (index_t d = 0; d < dim; d++)
    Y[k*dim+d] = topology.points()[k][d];
  domain_searcher->set_points( topology.points().nb() , Y.data() );

  printf("computing cells\n");
  for (index_t k = 0; k < nb_points; k++) {

    // find the closest vertex in the domain to this site
    real_t distance;
    index_t nearest;
    domain_searcher->get_nearest_neighbors( 1 , sites[k] , &nearest , &distance );
    //printf("nearest = %lu at distance %g\n",nearest,distance);

    std::vector<index_t> ball;
    topology.inverse().ball( nearest , ball );
    //print_inline(ball,"ball");

    voronoi::Cell cell( k , sites , topology , *searcher );
    cell.compute(ball);
    diagram.add_cell( cell );
  }
  printf("done computing cells\n");
  diagram.create_field();

  printf("volume = %g\n",diagram.volume());

  graphics::Viewer vis;
  vis.add(diagram);

  vis.run(AVRO_FULL_UNIT_TEST);
}
UT_TEST_CASE_END( test_2d )


UT_TEST_SUITE_END( voronoi_cell_test_suite )
