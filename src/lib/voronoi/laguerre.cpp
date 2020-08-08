#include "voronoi/laguerre.h"

#include "voronoi/voronoi_cell.h"

namespace avro
{

namespace delaunay
{

LaguerreDiagram::LaguerreDiagram( const Points& sites , const Topology<Simplex>& domain , const std::vector<real_t>& weights ) :
  Topology<Polytope>(points_,domain.number()),
  points_(sites.dim()),
  domain_points_(sites.dim()+1),
  domain_(domain_points_,domain.number()),
  weight_(weights),
  delaunay_(sites.dim()+1),
  exact_(true)
{
  avro_assert( sites.nb() > 0 );
  coord_t dim = sites.dim();

  // copy the domain, considering that it is embedded into dim + 1
  for (index_t k=0;k<domain.points().nb();k++)
  {
    domain_points_.create( domain.points()[k] );
    domain_points_[k][dim] = 0.0;
  }
  domain_.TopologyBase::copy( domain );

  // compute the maximum weight
  if (weight_.size() == 0)
    weight_.resize( sites.nb() , 0. );
  real_t wm = * std::max_element( weight_.begin() , weight_.end() ) + 1e-16;

  for (index_t k=0;k<sites.nb();k++)
  {
    delaunay_.create( sites[k] );
    delaunay_[k][dim] = std::sqrt( wm - weight_[k] );
  }
}

void
LaguerreDiagram::compute()
{
  // compute the voronoi diagram
  VoronoiDiagram diagram( delaunay_ , domain_ , true );
  diagram.compute(exact_);

  // copy the points, removing the extra coordinate (in dim+1)
  for (index_t k=0;k<diagram.points().nb();k++)
  {
    points_.create( diagram.points()[k] );
    std::vector<int> bisectors = diagram.points().incidence().get(k);
    points_.incidence().add( bisectors.data() , bisectors.size() );
  }
  TopologyBase::copy( diagram );
}

} // delaunay

} // avro
