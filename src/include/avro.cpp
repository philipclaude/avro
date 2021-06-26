#include "avro.h"

#include "adaptation/adapt.h"
#include "adaptation/metric.h"
#include "adaptation/parallel.h"

#include "element/simplex.h"

#include "geometry/entity.h"
#include "geometry/model.h"

#include "library/ckf.h"
#include "library/factory.h"

#include "mesh/facets.h"
#include "mesh/mesh.h"
#include "mesh/partition.h"
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
{}

Context::Context( const Context& ctx ) :
  number_(ctx.number()),
  dim_(ctx.dim()),
  udim_(ctx.udim()),
  model_(ctx.model()),
  points_(nullptr),
  topology_(nullptr)
{
  //parameters_.curved() = ctx.parameters().curved();
  import_model();
}

index_t
Context::nb_bodies() const
{
  avro_assert( model_ != nullptr );
  return model_->nb_bodies();
}

int
Context::facet_geometry( const index_t* v , index_t nv ) const
{
  Entity* entity = BoundaryUtils::geometryFacet( *points_.get() , v , nv );
  avro_assert( entity != nullptr );
  return entity2id_.at(entity);
}

void
Context::define_geometry( const std::string& geometry )
{
  // load the model using the factory
  bool curved;
  model_ = library::get_geometry(geometry,curved);
  parameters_.set( "curved" , curved );

  // import the data from the model
  import_model();
}

void
Context::define_geometry( const Model& model )
{
  avro_implement;
}

void
Context::import_model()
{
  // save the entities
  std::vector<Entity*> entities;
  model_->get_entities(entities);

  index_t number = 0;
  for (index_t k=0;k<entities.size();k++)
    number = (number > entities[k]->number()) ? number : entities[k]->number();

  std::vector<index_t> count(number+1,0);
  for (index_t k=0;k<entities.size();k++)
    count[entities[k]->number()]++;
  //print_inline(count,"count");

  std::vector<index_t> offset(number+1,0);
  for (index_t k=1;k<=number;k++)
    offset[k] = offset[k-1] + count[k-1];
  //print_inline(offset,"offset");

  std::vector<index_t> ids(entities.size());
  std::vector<index_t> recount(number+1,0);
  for (index_t k=0;k<entities.size();k++)
  {
    ids[k] = offset[entities[k]->number()] + recount[ entities[k]->number() ];
    recount[entities[k]->number()]++;
  }
  //print_inline(recount,"recount");

  for (index_t k=0;k<entities.size();k++)
  {
    id2entity_.insert( {ids[k] , entities[k]} );
    entity2id_.insert( {entities[k],ids[k]} );
    //entities[k]->print_header();
    //printf("--> maps to id %lu\n",ids[k]);
  }
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
  load_coordinates(x);
  load_simplices(s);
}

void
Context::load_coordinates( const std::vector<real_t>& x )
{
  index_t nb_points = x.size() / dim_;
  if (points_ != nullptr)
    points_->clear();
  points_ = std::make_shared<Points>(dim_,udim_);

  for (index_t k=0;k<nb_points;k++)
    points_->create( &x[k*dim_] );
}

void
Context::load_simplices( const std::vector<index_t>& s )
{
  index_t nb_simplices = s.size() / index_t(number_+1);

  if (topology_ != nullptr)
    topology_->clear();
  topology_ = std::make_shared<Topology<Simplex>>(*points_.get(),number_);

  for (index_t k=0;k<nb_simplices;k++)
    topology_->add( &s[(number_+1)*k] , number_+1 );
}

void
Context::load_local2global( const std::vector<index_t>& local2global )
{
  avro_assert( points_ != nullptr );
  avro_assert( points_->nb() == local2global.size() );

  for (index_t k = 0; k < local2global.size(); k++)
    points_->set_global( k , local2global[k] );
}

Entity*
Context::id2geometry( int id ) const
{
  if (id < 0) return nullptr;
  if (id2entity_.find(id) == id2entity_.end())
  {
    for (std::map<int,Entity*>::const_iterator it=id2entity_.begin();it!=id2entity_.end();++it)
    {
      it->second->print_header();
      printf("--> with id = %d\n",it->first);
    }
    printf("could not find entity with id = %d\n",id);
  }
  avro_assert( id2entity_.find(id) != id2entity_.end() );
  return id2entity_.at(id);
}

void
Context::get_geometry_ids( std::vector<int>& ids ) const
{
  ids.clear();
  for (std::map<int,Entity*>::const_iterator it=id2entity_.begin();it!=id2entity_.end();++it)
    ids.push_back( it->first );
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
  std::vector<symd<real_t>> metric;
  index_t count = 0;
  for (index_t k=0;k<nb_metric;k++)
  {
    symd<real_t> mk(number_);
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
  AdaptationParameters params(parameters_);
  //params.prefix() = parameters_.get_param<std::string>("prefix");
  //params.write_conformity() = parameters_.get_param<bool>("write conformity");
  //params.write_mesh() = parameters_.get_param<bool>("write mesh");
  //params.output_redirect() = parameters_.get_param<std::string>("output redirect");

  AdaptationProblem problem = {mesh_in,metric,params,mesh_out};
  int result = ::avro::adapt<Simplex>( problem );
  if (result != 0) return 1;

  // copy the adapted mesh into the context data
  mesh_out.points().copy( *points_.get() );
  topology_->TopologyBase::copy( mesh_out.retrieve<Simplex>(0) );

  return 0; // success
}

int
Context::adapt_parallel( const std::vector<real_t>& m )
{
  avro_assert_msg( points_ != nullptr , "points are not defined" );
  avro_assert_msg( topology_ != nullptr , "topology is not defined" );

  index_t nb_rank = (number_+1)*number_/2;
  index_t nb_metric = m.size() / nb_rank;
  avro_assert( nb_metric == points_->nb() );

  // save the metrics in the proper format
  std::vector<symd<real_t>> metric;
  index_t count = 0;
  for (index_t k=0;k<nb_metric;k++)
  {
    symd<real_t> mk(number_);
    for (index_t i=0;i<number_;i++)
    for (index_t j=i;j<number_;j++)
      mk(i,j) = m[count++];
    metric.push_back(mk);
  }

  // setup the adaptation manager
  AdaptationParameters params(parameters_);
  params.set( "partitioned" , true );
  params.set( "allow serial" , true );
  params.set( "insertion volume factor" ,  -1.0 );
  params.set( "curved" , false);
  params.set( "limit metric" , true );
  params.set( "max parallel passes" , index_t(3) );
  params.set( "elems per processor" , index_t(5000) );
  params.set("has uv", true);
  params.set( "swapout" , false);

  std::vector<Entity*> entities;
  model_->get_tessellatable_entities(entities);

  topology_->build_structures();
  AdaptationManager<Simplex> manager( *topology_.get() , metric , params );

  manager.topology().set_entities( entities );
  printf("there are %lu entities\n",entities.size());

  manager.adapt();

  // copy the adapted mesh into the context data
  topology_->TopologyBase::copy( manager.topology() );
  manager.topology().points().copy( *points_.get() );

  mpi::barrier();

  return 0; // success
}

void
Context::partition() {

  AdaptationParameters params;
  params.set( "partitioned" , false );

  // dummy metrics to initialize the adaptation manager
  // yes, this means there will be an unncessary synchronization of metric data
  // which will then get thrown away, but it should be a low overhead since this
  // should only be done **once** at the beginning of a simulation
  index_t nb_points = points_->nb();
  std::vector<symd<real_t>> metric( nb_points );

  topology_->build_structures();
  AdaptationManager<Simplex> manager( *topology_.get() , metric , params );

  // copy the partitioned mesh into the context data
  topology_->TopologyBase::copy( manager.topology() );
  manager.topology().points().copy( *points_.get() );
}

void
Context::retrieve_mesh( std::vector<real_t>& x , std::vector<index_t>& s ) const
{
  avro_assert_msg( points_ != nullptr , "points are not defined" );
  avro_assert_msg( topology_ != nullptr , "topology is not defined" );

  x.resize( points_->nb()*points_->dim() );
  index_t i = 0;
  for (index_t k=0;k<points_->nb();k++)
  for (coord_t d=0;d<points_->dim();d++)
    x[i++] = (*points_)(k,d);

  index_t nv = number_+1;
  s.resize( topology_->nb()*nv );
  i = 0;
  for (index_t k=0;k<topology_->nb();k++)
  for (index_t d=0;d<nv;d++)
    s[i++] = (*topology_)(k,d);
}

void
Context::retrieve_local2global( std::vector<index_t>& local2global ) const
{
  avro_assert_msg( points_ != nullptr , "points are not defined" );

  local2global.resize( points_->nb() );
  for (index_t k = 0; k < points_->nb(); k++)
    local2global[k] = points_->global(k);
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
    g[k] = entity2id_.at(points_->entity(k));
    for (coord_t d=0;d<udim_;d++)
      u[k*udim_+d] = points_->u(k,d);
  }
}

void
Context::retrieve_boundary( std::vector<std::vector<index_t>>& facets ,
                            std::vector<int>& geometry , bool interior ) const
{
  avro_assert_msg( topology_ != nullptr , "topology is not defined" );
  Boundary<Simplex> boundary(*topology_.get());
  boundary.extract(interior);

  facets.resize( boundary.nb_children() , std::vector<index_t>() );
  geometry.resize( boundary.nb_children() );
  for (index_t k=0;k<boundary.nb_children();k++)
  {
    Entity* entity = boundary.entity(k);
    const Topology<Simplex>& bk = boundary.child(k);
    avro_assert( entity->number() == index_t(number_-1) );
    avro_assert( bk.number() == entity->number() );
    for (index_t j=0;j<bk.nb();j++)
    for (index_t i=0;i<bk.nv(j);i++)
      facets[k].push_back( bk(j,i) );
    geometry[k] = entity2id_.at(entity);
  }
}

void
Context::retrieve_boundary_parallel( std::vector<std::vector<index_t>>& faces ,
                                     std::vector<int>& geometry ) const
{
  avro_assert_msg( topology_ != nullptr , "topology is not defined" );

  std::map<Entity*,int>::const_iterator it;
  std::map<Entity*,int> entity2index;
  for (it = entity2id_.begin(); it != entity2id_.end(); ++it) {
    if (it->first->number() == number_-1) {
      entity2index.insert( {it->first,entity2index.size()} );
    }
  }
  index_t nb_boundary = entity2index.size();

  faces.resize( nb_boundary );
  geometry.resize( nb_boundary );

  // compute the boundary
  Facets facets(*topology_.get());
  facets.compute();
  std::vector<index_t> facet(topology_->number());
  for (index_t k = 0; k < facets.nb(); k++) {

    // skip partition interface facets or interior facets
    // (that do not lie on wakes, i.e. interior geometries)
    facets.retrieve(k,facet);
    Entity* entity = BoundaryUtils::geometryFacet( topology_->points() , facet.data() , facet.size());
    if (entity == nullptr) continue;

    // this is a geometry facet that is on either an external or interior geometry
    index_t idx = entity2index.at(entity);
    for (index_t j = 0; j < facet.size(); j++)
      faces[idx].push_back(facet[j]);
  }

  for (it = entity2index.begin(); it != entity2index.end(); ++it)
    geometry[it->second] = entity2id_.at(it->first);
}

}
