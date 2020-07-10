#if 0
#include "adaptation/parallel0.hpp"
#else
#include "common/mpi.hpp"
#include "common/process.h"

#include "adaptation/adapt.h"
#include "adaptation/parallel.h"
#include "adaptation/parameters.h"

#include "element/simplex.h"

#include "library/meshb.h"

#include "mesh/boundary.h"
#include "mesh/facets.h"
#include "mesh/mpi_tags.h"
#include "mesh/partition.h"
#include "mesh/points.h"
#include "mesh/topology.h"

#include "numerics/matrix.h"

namespace avro
{

#ifdef AVRO_MPI

#define TAG_METRIC_INDICES 301
#define TAG_METRIC_ENTRIES 302

#define TAG_FIXED_GLOBALS 401
#define TAG_GLOBAL_KEYS 402
#define TAG_GLOBAL_VALUES 403
#define TAG_POINT_OFFSET 404

void
pprintf(const char* format , ... )
{
  va_list args;
  va_start(args,format);
  index_t rank = mpi::rank();
  printf("[processor %2lu]: ",rank);
  printf(format,args);
  va_end(args);
}

template<typename type>
AdaptationManager<type>::AdaptationManager( const Topology<type>& topology , const std::vector<VertexMetric>& metrics , AdaptationParameters& params ) :
  params_(params),
  topology_( topology.points().dim() , topology.points().udim() , topology.number() )
{
  // save some mpi stuff
  rank_ = mpi::rank();

  // initialize our working topology (may require partitioning)
  initialize(topology,metrics);
}

template<typename type>
void
AdaptationManager<type>::send_metrics( index_t receiver , const std::vector<index_t>& indices , const std::vector<VertexMetric>& metrics , bool global_flag )
{
  // number of unique entries in a metric
  index_t nb_entry = topology_.number()*(topology_.number()+1)/2;
  std::vector<index_t> global_indices( indices.size() , 0 );
  std::vector<real_t> metric_data( indices.size()*nb_entry , 0.0 );

  index_t entry = 0;
  if (global_flag)
  {
    // indices are global, we we need to retrieve the local metrics
    for (index_t k=0;k<indices.size();k++)
    {
      index_t g = indices[k];
      const VertexMetric& m = metrics[g];
      for (index_t i=0;i<m.n();i++)
      for (index_t j=i;j<m.n();j++)
        metric_data[entry++] = m(i,j);
    }
    global_indices = indices;
  }
  else
  {
    // indices are local, so we need to convert to global
    for (index_t k=0;k<indices.size();k++)
    {
      const VertexMetric& m = metrics[indices[k]];
      for (index_t i=0;i<m.n();i++)
      for (index_t j=i;j<m.n();j++)
        metric_data[entry++] = m(i,j);
      global_indices[k] = topology_.points().global( indices[k] );
    }
  }

  // send the data
  mpi::send( mpi::blocking{} , global_indices , receiver , TAG_METRIC_INDICES );
  mpi::send( mpi::blocking{} , metric_data , receiver , TAG_METRIC_ENTRIES );
}

template<typename type>
void
AdaptationManager<type>::receive_metrics( index_t sender , bool overwrite )
{
  // number of unique entries in a metric
  index_t n = topology_.number()*(topology_.number()+1)/2;

  // option to overwrite means all metrics will be cleared
  if (overwrite)
    metrics_.resize( topology_.points().nb() , VertexMetric(topology_.number()) );

  // receive the data
  std::vector<index_t> global_indices = mpi::receive< std::vector<index_t> >( sender , TAG_METRIC_INDICES );
  std::vector<real_t> metric_data = mpi::receive< std::vector<real_t> >( sender , TAG_METRIC_ENTRIES );

  // convert the global_indices to a map so we can look things up faster
  std::map<index_t,index_t> global;
  for (index_t k=0;k<global_indices.size();k++)
    global.insert( {global_indices[k] , k} );

  for (index_t k=0;k<topology_.points().nb();k++)
  {
    std::map<index_t,index_t>::iterator it = global.find( topology_.points().global(k) );
    avro_assert( it!=global.end() );

    // determine which index this is in the list
    index_t idx = it->second;

    // unpack the metric
    VertexMetric m( topology_.number() );
    index_t entry = 0;
    for (index_t i=0;i<m.n();i++)
    for (index_t j=i;j<m.n();j++)
      m(i,j) = metric_data[ idx*n + entry++ ];
    metrics_[k] = m;
  }
}

template<typename type>
void
AdaptationManager<type>::initialize(const Topology<type>& topology , const std::vector<VertexMetric>& metrics )
{
  index_t nb_rank = mpi::size();
  index_t nb_partition = nb_rank;

  if (!params_.partitioned())
  {
    // extract the entities so we can assign them to the partitions
    std::vector<Entity*> entities;
    for (index_t k=0;k<topology.points().nb();k++)
    {
      if (topology.points().entity(k)==nullptr) continue;
      entities.push_back( topology.points().entity(k) );
    }
    uniquify(entities);

    if (rank_ == 0)
    {
      // partition the input topology
      Partition<type> partition(topology);
      partition.compute( nb_partition );

      // assume the mesh is balanced now
      params_.balanced() = true;

      // extract the partitions
      std::vector<std::shared_ptr<Topology_Partition<type>>> pieces(nb_partition);
      partition.get(pieces);

      // steal partition 0 with the metrics into our working partition
      topology_.TopologyBase::copy( *pieces[0] );
      pieces[0]->points().copy( topology_.points() );
      pieces[0]->set_entities( entities );
      metrics_.resize( pieces[0]->points().nb() );
      for (index_t j=0;j<pieces[0]->points().nb();j++)
        metrics_[j] = metrics[ pieces[0]->points().global(j) ];

      // send off the remaining partitions to the processors
      for (index_t k=1;k<nb_rank;k++)
      {
        pieces[k]->set_entities(entities);
        pieces[k]->send(k);

        std::vector<index_t> indices(pieces[k]->points().nb());
        for (index_t j=0;j<pieces[k]->points().nb();j++)
          indices[j] = pieces[k]->points().global(j);
        send_metrics( k , indices , metrics , true );
      }
    }
    else
    {
      // receive our partition
      pprintf("waiting to receive mesh\n");
      topology_.set_entities(entities);
      topology_.receive(0);

      // receive our set of metrics
      receive_metrics(0,true);
    }
  }
  else
  {
    // copy the input topology into the working topology
    topology_.TopologyBase::copy( topology );
    topology.points().copy( topology_.points() );

    // copy the input metrics
    metrics_ = metrics;
  }
  mpi::barrier();

  // balance the partitions so they have an equal adaptation work
  if (!params_.balanced())
  {
    // everything should still be saved into the partition
    balance();
  }
}

template<typename type>
void
AdaptationManager<type>::balance(real_t alpha , real_t beta)
{
  if (alpha<0) alpha = 1.0;
  if (beta<0) beta = 1.0;

  //avro_assert( metrics_.size() == topology_.points().nb() );

  // compute the weights on all the elements in our partition

  // compute the adjacency information

  // call the repartitioner/load-balancer

  // transfer data to another partition
  // TODO topology and metrics

  // receive data from another partition
  // TODO topology and metrics

  mpi::barrier();
}

template<typename type>
void
AdaptationManager<type>::migrate_parmetis()
{
  index_t nb_rank = mpi::size();

  // migrate the interface into the interior of the partitions
  // simultaneously performing a load balance
  topology_.build_structures();

  // assign global identifiers for all elements

  // exchange partition boundaries to all processors
  PartitionBoundary<type> boundary(topology_.number()-1);
  boundary.compute(topology_,rank_);

  PartitionBoundary<type> interface( topology_.number()-1 );
  if (rank_ == 0)
  {
    // receive and accumulate all boundary information
    interface.append( boundary );
    for (index_t k=1;k<nb_rank;k++)
    {
      PartitionBoundary<type> boundary_k( topology_.number()-1 );
      boundary_k.receive(k);
      interface.append(boundary_k);

      printf("nb_facets = %lu\n",boundary_k.nb());
    }
    printf("total number of interface facets = %lu\n",interface.nb());
  }
  else
  {
    // send our partition boundary information
    printf("sending boundary with %lu facets\n",boundary.nb());
    boundary.send(0);
  }
  mpi::barrier();

  if (rank_ == 0 )
  {
    // send the accumulated boundary information
    printf("sending interface with %lu facets\n",interface.nb());
    for (index_t k=1;k<nb_rank;k++)
      interface.send(k);
  }
  else
  {
    // receive the accumulated boundary information
    interface.receive(0);
  }
  mpi::barrier();

  // fill in the missing elemR/partR information
  //boundary.fill( interface );

  // build up the local adjacency graph, accounting for elements (graph vertices)
  // on other processors stored in the elemR/partR information

  // determine a weight on each element using the metric

  // build up the graph representation for this partition
  // and weights the vertices of the graph (elements of the mesh)

  // call parmetis to perform the load balance

  // retrieve which elements we keep and which we send away
  // also send away the metrics

  // receive elements from another partition
  // also receive the metrics
}

template<typename type>
void
AdaptationManager<type>::migrate_native()
{
  // setup a buddy system
  index_t buddy = 0;
  if (rank_ == 0)
  {
    // build up the interprocessor graph
    // with vertex weights defined be current processor work

    // compute a matching of this graph

    // send off the processor pairs
  }
  else
  {
    // wait to receive our processor pair
  }
  mpi::barrier();

  pprintf("my buddy is %lu\n",buddy);
}

template<typename type>
void
AdaptationManager<type>::migrate()
{
  // TODO add parameter to switch between parmetis/native
  migrate_parmetis();
}

void
get_fixed_globals( const Points& points , std::vector<index_t>& global )
{
  for (index_t k=0;k<points.nb();k++)
  {
    if (!points.fixed(k)) continue;
    global.push_back(points.global(k));
  }
}

template<typename key,typename value>
struct map_builder : public std::map<key,value>
{
  map_builder( const std::vector<key>& keys , const std::vector<value>& values )
  {
    avro_assert( keys.size() == values.size() );
    for (index_t k=0;k<keys.size();k++)
      this->insert( {keys[k] , values[k] } );
  }
};


template<typename type>
void
AdaptationManager<type>::synchronize()
{
  index_t nb_rank = mpi::size();

  // synchronize all the global point indices between processors
  // the root processor needs to know how many vertices there are in total
  std::vector<index_t> global;
  get_fixed_globals( topology_.points() , global );
  std::map<index_t,index_t> fixed_map;
  index_t nb_points_total = 0;
  std::vector<index_t> nb_points( nb_rank , 0 );
  std::vector<index_t> pt_offset( nb_rank , 0 );
  if (rank_ == 0)
  {
    for (index_t j=0;j<global.size();j++)
      fixed_map.insert( {global[j] , j } );

    // receive the total number of vertices from all processors
    nb_points[0] = topology_.points().nb();

    nb_points_total = nb_points[0];
    std::vector<index_t> nb_interior( nb_rank , 0 );
    nb_interior[0] = nb_points[0] - global.size();
    for (index_t k=1;k<nb_rank;k++)
    {
      // receive the total number of points
      nb_points[k] += index_t(mpi::receive<int>( k , TAG_MISC ));
      nb_points_total += nb_points[k];

      // receive the list of global indices
      std::vector<index_t> globals_k = mpi::receive<std::vector<index_t>>( k , TAG_FIXED_GLOBALS );
      nb_interior[k] = nb_points[k] - globals_k.size();

      for (index_t j=0;j<globals_k.size();j++)
      {
        index_t g = globals_k[j];
        if (fixed_map.find(g)!=fixed_map.end()) continue;
        fixed_map.insert( {g,fixed_map.size()} );
      }
    }

    // compute the point offset for every processor
    pt_offset[0] = fixed_map.size();
    for (index_t k=1;k<nb_rank;k++)
    {
      pt_offset[k] = pt_offset[k-1] + nb_interior[k-1];
    }

  }
  else
  {
    // send the number of points on this processor
    mpi::send( mpi::blocking{} , int(topology_.points().nb()) , 0 , TAG_MISC );

    // send list of global indices
    mpi::send( mpi::blocking{} , global , 0 , TAG_FIXED_GLOBALS );
  }
  mpi::barrier();

  if (rank_ == 0)
  {
    // assign the new global indices for this processor
    index_t count = fixed_map.size();
    for (index_t j=0;j<topology_.points().nb();j++)
    {
      // fixed points should still be at the beginning!
      if (!topology_.points().fixed(j))
      {
        topology_.points().set_global( j , count++ );
        continue;
      }

      index_t g = topology_.points().global(j);
      topology_.points().set_global( j , fixed_map[g] );
    }

    // send global indices to all processors
    std::vector<index_t> keys,values;
    for (std::map<index_t,index_t>::const_iterator it=fixed_map.begin();it!=fixed_map.end();++it)
    {
      keys.push_back( it->first );
      values.push_back( it->second );
    }
    for (index_t k=1;k<nb_rank;k++)
    {
      mpi::send( mpi::blocking{} , keys , k , TAG_GLOBAL_KEYS );
      mpi::send( mpi::blocking{} , values , k , TAG_GLOBAL_VALUES );
      mpi::send( mpi::blocking{} , pt_offset[k] , k , TAG_POINT_OFFSET );
    }
  }
  else
  {
    // receive the old global indices
    std::vector<index_t> keys = mpi::receive<std::vector<index_t>>(0,TAG_GLOBAL_KEYS);

    // receive new global indices
    std::vector<index_t> values = mpi::receive<std::vector<index_t>>(0,TAG_GLOBAL_VALUES);

    // receive the point offset
    index_t offset = mpi::receive<index_t>(0,TAG_POINT_OFFSET);

    // make a map!
    map_builder<index_t,index_t> global_map(keys,values);

    index_t count = 0;
    for (index_t k=0;k<topology_.points().nb();k++)
    {
      // fixed points should still be at the beginning!
      if (!topology_.points().fixed(k))
      {
        topology_.points().set_global( k , offset + count++ );
        continue;
      }

      // current global index
      index_t g = topology_.points().global(k);
      avro_assert( global_map.find(g)!=global_map.end() );

      // assign the new global index
      topology_.points().set_global( k , global_map.at(g) );
    }
  }

  mpi::barrier();
}

template<typename type>
bool
AdaptationManager<type>::analyze()
{
  bool result = false;

  // compute metric conformity on our processor

  // communicate the result to all processors

  // if a single processor is unhappy, then we are not done

  mpi::barrier();

  return result;
}

template<typename type>
void
AdaptationManager<type>::fix_boundary()
{
  // first unfix all the vertices
  for (index_t k=0;k<topology_.points().nb();k++)
    topology_.points().set_fixed(k,false);

  // compute the boundary
  Facets facets(topology_);
  facets.compute();
  std::vector<index_t> facet(topology_.number());
  std::vector<index_t> pts;
  index_t nb_fixed_facets = 0;
  index_t nb_bnd = 0;
  for (index_t k=0;k<facets.nb();k++)
  {
    if (!facets.boundary(k)) continue;
    nb_bnd++;
    facets.retrieve(k,facet);
    Entity* entity = BoundaryUtils::geometryFacet( topology_.points() , facet.data() , facet.size());
    if (entity!=nullptr) continue;
    nb_fixed_facets++;

    for (index_t j=0;j<facet.size();j++)
      pts.push_back(facet[j]);
  }
  uniquify(pts);
  std::sort(pts.begin(),pts.end());

  if (mpi::size() == 1 )
    avro_assert( nb_fixed_facets == 0 );

  // fix all the partition boundary points
  for (index_t j=0;j<pts.size();j++)
    topology_.points().set_fixed(pts[j],true);

  // move all the fixed points to the beginning of the list
  // this should also adjust the global indices
  // but we'll need to adjust the metrics
  std::map<index_t,index_t> point_map;
  topology_.move_to_front( pts , &point_map );

  // map the metrics
  coord_t number = topology_.number();
  std::vector<VertexMetric> mapped_metrics(metrics_.size(),VertexMetric(number));
  for (index_t k=0;k<metrics_.size();k++)
    mapped_metrics[ point_map[k] ] = metrics_[k];
  metrics_ = mapped_metrics;
}

template<typename type>
void
AdaptationManager<type>::adapt()
{
  mpi::communicator& comm = ProcessMPI::get_comm();
  UNUSED(comm);
  coord_t number = topology_.number();

  bool done = false;
  for (index_t pass=0;pass<params_.max_passes();pass++)
  {
    // fix the boundary of the topology
    // and move the fixed points to the beginning of the points structure
    fix_boundary();

    // do the adaptation
    // create the mesh we will write to
    // create a mesh and add the topology
    Mesh mesh(number,number);
    std::shared_ptr<Topology<type>> ptopology = std::make_shared<Topology<type>>(mesh.points(),number);
    ptopology->TopologyBase::copy(topology_);
    mesh.add(ptopology);
    topology_.points().copy(mesh.points());
    Mesh mesh_out(number,number);

    //if (rank_ > 0)
      params_.output_redirect() = "adaptation-output-proc"+std::to_string(rank_)+".txt";
    params_.export_boundary() = false;
    params_.prefix() = "mesh-proc"+std::to_string(rank_);
    AdaptationProblem problem = {mesh,metrics_,params_,mesh_out};

    try
    {
      #if 1
      pprintf("adapting mesh\n");
      ::avro::adapt<type>( problem );

      // clear the topology and copy in the output topology
      topology_.clear();
      topology_.points().clear();
      topology_.TopologyBase::copy( mesh_out.template retrieve<type>(0) );
      mesh_out.points().copy( topology_.points() );
      #endif

      // unset any global indices in unfixed points, since these need to be renumbered
      for (index_t k=0;k<topology_.points().nb();k++)
      {
        if (topology_.points().fixed(k)) continue;
        topology_.points().set_global( k , 0 );
      }
    }
    catch(...)
    {
      mpi::abort(1);
    }
    mpi::barrier();

    // synchronize all the global point indices
    synchronize();

    // analyze whether we are done
    done = analyze();
    if (done) break;

    // migrate the interface between partitions
    migrate();

    break; // only one pass for now..
  }

  // balance the mesh between the partitions
  balance();
}

template<typename type>
class PartitionField : public Field<type,real_t>
{
public:
  PartitionField( const Topology<type>& topology ) :
    Field<type,real_t>(topology,1,DISCONTINUOUS)
  {
    this->element().set_basis( BasisFunctionCategory_Lagrange );
  }

  index_t nb_rank() const override { return 1; }
  std::string get_name( index_t j ) const override
  {
    return std::to_string(j);
  }
};

template<typename type>
index_t
AdaptationManager<type>::get_total_nb_points() const
{
  index_t nb_rank = mpi::size();

  // get the maximum global index
  index_t nb_points = 0;
  for (index_t k=0;k<topology_.points().nb();k++)
  {
    if (topology_.points().global(k) > nb_points)
      nb_points = topology_.points().global(k);
  }

  if (rank_ == 0)
  {
    for (index_t k=1;k<nb_rank;k++)
    {
      index_t nb_points_k = index_t(mpi::receive<int>( k , TAG_MISC ));
      nb_points = std::max(nb_points,nb_points_k);
    }
  }
  else
  {
    mpi::send( mpi::blocking{} , int(nb_points) , 0 , TAG_MISC );
  }
  mpi::barrier();
  return nb_points;
}

template<typename type>
void
AdaptationManager<type>::append_partition( index_t p , Topology<type>& topology , Topology_Partition<type>& partition ,
                                           std::vector<bool>& created , std::vector<real_t>& partition_index ) const
{
  for (index_t j=0;j<partition.points().nb();j++)
  {
    // retrieve the global index of this point
    index_t g = partition.points().global(j);
    if (created[g]) continue;
    created[g] = true;
    for (index_t d=0;d<topology.points().dim();d++)
      topology.points()[g][d] = partition.points()[j][d];
  }

  for (index_t j=0;j<partition.nb();j++)
  {
    partition_index.push_back(p);
    std::vector<index_t> s = partition.get(j);
    for (index_t i=0;i<s.size();i++)
    {
      s[i] = partition.points().global(s[i]);
      avro_assert( s[i] < topology.points().nb() );
    }
    topology.add(s.data(),s.size());
  }
}

template<typename type>
void
AdaptationManager<type>::retrieve( Topology<type>& topology )
{
  // fill the topology with the elements
  // fill the points with the vertices
  // assign the globals into the points
  index_t nb_rank = mpi::size();
  index_t nb_points = get_total_nb_points() + 1; // +1 because of 0-bias

  if (rank_==0)
  {
    std::vector<bool> created(nb_points,false);

    // initialize all the points (data filled below)
    std::vector<real_t> zero(topology.points().dim());
    for (index_t j=0;j<nb_points;j++)
      topology.points().create(zero.data());

    // append the root partition into the output topology
    std::vector<real_t> partition_index;
    append_partition( 0 , topology , topology_ , created , partition_index );
    for (index_t k=1;k<nb_rank;k++)
    {
      // receive the partition from processor k and append it into the output topology
      Topology_Partition<type> partition( topology_.points().dim() , topology_.points().udim() , topology_.number() );
      partition.receive(k);
      append_partition( k , topology , partition , created , partition_index );
    }

    // fill the partition field
    field_ = std::make_shared<PartitionField<type>>(topology);
    topology.element().set_basis( BasisFunctionCategory_Lagrange );
    field_->build();
    for (index_t k=0;k<topology.nb();k++)
    for (index_t j=0;j<topology.nv(k);j++)
      (*field_)(k,j) = partition_index[k];
    topology.fields().make("partition",field_);

  }
  else
  {
    // send the working topology to the root processor
    topology_.send(0);
  }
  mpi::barrier();
}

template class AdaptationManager<Simplex>;

#endif

} // avro

#endif
