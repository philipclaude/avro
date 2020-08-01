#include "avro.h"

#include "adaptation/adapt.h"
#include "adaptation/metric.h"
#include "adaptation/parallel.h"

#include "element/simplex.h"

#include "geometry/entity.h"
#include "geometry/model.h"

#include "library/ckf.h"
#include "library/factory.h"

#include "mesh/mesh.h"
#include "mesh/points.h"
#include "mesh/topology.h"

namespace avro
{

Context::Context( coord_t number , coord_t dim , coord_t udim ) :
  number_(number),
  dim_(dim),
  udim_(udim),
  model_(nullptr),
  points_(nullptr),
  topology_(nullptr)
{
  parameters_.standard();
}

void
Context::define_geometry( const std::string& geometry )
{
  // load the model using the factory
  bool curved;
  model_ = library::get_geometry(geometry,curved);
  parameters_.curved() = curved;

  // save the entities
  model_->get_entities(entities_);
}

void
Context::define_mesh( const std::string& mesh )
{
  coord_t number = number_;
  std::shared_ptr<Mesh> mesh_ptr;
  std::shared_ptr<TopologyBase> topology;
  mesh_ptr = library::get_mesh( mesh , topology , number );
  avro_assert_msg( number == number_ , "mismatch in requested mesh topological number" );
  load_mesh( topology->points().data() , topology->data() );
}

void
Context::attach_geometry()
{
  avro_assert( points_ != nullptr );
  avro_assert( model_ != nullptr );
  points_->attach(*model_.get());
}

void
Context::load_mesh( const std::vector<real_t>& x , const std::vector<index_t>& s )
{
  index_t nb_points = x.size() / dim_;
  index_t nb_simplices = s.size() / index_t(number_+1);

  if (points_ != nullptr)
    points_->clear();
  points_ = std::make_shared<Points>(dim_,udim_);

  for (index_t k=0;k<nb_points;k++)
    points_->create( &x[k*dim_] );

  if (topology_ != nullptr)
    topology_->clear();
  topology_ = std::make_shared<Topology<Simplex>>(*points_.get(),number_);

  for (index_t k=0;k<nb_simplices;k++)
    topology_->add( &s[(number_+1)*k] , number_+1 );
}

Entity*
Context::id2geometry( int id ) const
{
  for (index_t k=0;k<entities_.size();k++)
  {
    if (entities_[k]->identifier() == id)
      return entities_[k];
  }
  return nullptr;
}

void
Context::load_geometry( const std::vector<int>& g , const std::vector<real_t>& u )
{
  avro_assert( points_ != nullptr );
  avro_assert( g.size() == points_->nb() );
  avro_assert( u.size() == udim_*points_->nb() );

  for (index_t k=0;k<g.size();k++)
  {
    points_->set_entity( k , id2geometry(g[k]) );
    points_->set_param( k , &u[k*udim_] );
  }
}

int
Context::adapt( const std::vector<real_t>& m )
{
  avro_assert_msg( points_ != nullptr , "points are not defined" );
  avro_assert_msg( topology_ != nullptr , "topology is not defined" );

  index_t nb_rank = (number_+1)*number_/2;
  index_t nb_metric = m.size() / nb_rank;
  avro_assert( nb_metric == points_->nb() );

  // save the metrics in the proper format
  std::vector<numerics::SymMatrixD<real_t>> metric;
  index_t count = 0;
  for (index_t k=0;k<nb_metric;k++)
  {
    numerics::SymMatrixD<real_t> mk(number_);
    for (index_t i=0;i<number_;i++)
    for (index_t j=i;j<number_;j++)
      mk(i,j) = m[count++];
    metric.push_back(mk);
  }

  Mesh mesh_in(dim_,number_);
  Mesh mesh_out(dim_,number_);

  points_->copy( mesh_in.points() );
  mesh_in.add( topology_ );

  // setup the problem
  AdaptationProblem problem = {mesh_in,metric,parameters_,mesh_out};
  int result = ::avro::adapt<Simplex>( problem );
  if (result != 0) return 1;

  // copy the adapted mesh into the context data
  mesh_out.points().copy( *points_.get() );
  topology_->TopologyBase::copy( mesh_out.retrieve<Simplex>(0) );

  return 0; // success
}

void
Context::retrieve_mesh( std::vector<real_t>& x , std::vector<index_t>& s ) const
{
  avro_assert_msg( points_ != nullptr , "points are not defined" );
  avro_assert_msg( topology_ != nullptr , "topology is not defined" );

  x.assign( points_->data().begin() , points_->data().end() );
  s.assign( topology_->data().begin() , topology_->data().end() );
}

void
Context::retrieve_geometry( std::vector<int>& g , std::vector<real_t>& u ) const
{
  avro_assert( points_ != nullptr );

  g.resize( points_->nb() , -1 );
  u.resize( udim_*points_->nb() , unset_value );

  for (index_t k=0;k<points_->nb();k++)
  {
    if (points_->entity(k) == nullptr) continue;
    g[k] = points_->entity(k)->identifier();
    for (coord_t d=0;d<udim_;d++)
      u[k*udim_+d] = points_->u(k,d);
  }
}

void
Context::retrieve_boundary( std::vector<std::vector<index_t>>& facets ,
                            std::vector<int>& geometry , bool interior ) const
{
  Boundary<Simplex> boundary(*topology_.get());
  boundary.extract(interior);

  facets.resize( boundary.nb_children() , std::vector<index_t>() );
  geometry.resize( boundary.nb_children() );
  for (index_t k=0;k<boundary.nb_children();k++)
  {
    Entity* entity = boundary.entity(k);
    const Topology<Simplex> bk = boundary.child(k);
    facets[k].assign( bk.data().begin() , bk.data().end() );
    geometry[k] = entity->identifier();
  }
}

}
