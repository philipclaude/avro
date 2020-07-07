#if 0
#include "adaptation/parallel0.hpp"
#else
#include "common/mpi.hpp"
#include "common/process.h"

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

template<typename type>
AdaptationManager<type>::AdaptationManager( const Topology<type>& topology , const std::vector<VertexMetric>& metrics , AdaptationParameters& params ) :
  params_(params),
  topology_( topology.points().dim() , topology.points().udim() , topology.number() )
{
  initialize(topology,metrics);

  // save some mpi stuff
  rank_ = mpi::rank();
}

template<typename type>
void
AdaptationManager<type>::send_metrics( index_t receiver , const std::vector<index_t>& indices , const std::vector<VertexMetric>& metrics , bool global_flag )
{
  // number of unique entries in a metric
  index_t N = topology_.number()*(topology_.number()+1)/2;
  std::vector<index_t> global_indices( indices.size() , 0 );
  std::vector<real_t> metric_data( indices.size()*N , 0.0 );

  index_t entry = 0;
  if (global_flag)
  {
    // indices are global, we we need to retrieve the local metrics
    std::set<index_t> global(indices.begin(),indices.end());
    for (index_t k=0;k<topology_.points().nb();k++)
    {
      index_t g = topology_.points().global(k);
      if (global.find(g)==global.end()) continue;

      const VertexMetric& m = metrics_[k];
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
      const VertexMetric& m = metrics_[indices[k]];
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
  index_t N = topology_.number()*(topology_.number()+1)/2;

  // option to overwrite means all metrics will be cleared
  if (overwrite)
    metrics_.resize( topology_.points().nb() );

  // receive the data
  std::vector<index_t> global_indices = mpi::receive< std::vector<index_t> >( sender , TAG_METRIC_INDICES );
  std::vector<real_t> metric_data = mpi::receive< std::vector<real_t> >( sender , TAG_METRIC_ENTRIES );

  // convert the global_indices to a set so we can look things up faster
  std::map<index_t,index_t> global;
  for (index_t k=0;k<global_indices.size();k++)
    global.insert( {global_indices[k] , k} );

  for (index_t k=0;k<topology_.points().nb();k++)
  {
    std::map<index_t,index_t>::iterator it = global.find( topology_.points().global(k) );
    if (it==global.end()) continue; // not setting this metric

    // determine which index this is in the list
    index_t idx = it->second;

    // unpack the metric
    VertexMetric m( topology_.number() );
    index_t entry = 0;
    for (index_t i=0;i<m.n();i++)
    for (index_t j=i;j<m.n();j++)
      m(i,j) = metric_data[ idx*N + entry++ ];
    metrics_[k] = m;
  }
}

template<typename type>
void
AdaptationManager<type>::initialize( const Topology<type>& topology , const std::vector<VertexMetric>& metrics )
{
  index_t nb_rank = mpi::size();

  if (!params_.partitioned())
  {
    if (rank_ == 0)
    {
      // partition the input topology
      Partition<type> partition(topology);
      partition.compute( nb_rank );

      // extract the partitions
      std::vector<std::shared_ptr<Topology_Partition<type>>> pieces;
      partition.get(pieces);

      // steal partition 0 into our working partition
      topology_.TopologyBase::copy( *pieces[0] );
      pieces[0]->points().copy( topology_.points() );

      // send off the remaining partitions to the processors
      for (index_t k=1;k<nb_rank;k++)
      {
        pieces[k]->send(k);

        std::vector<index_t> indices = linspace(pieces[k]->points().nb());
        send_metrics( k , indices , metrics , true );
      }
    }
    else
    {
      // receive our partition
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

  avro_assert( metrics_.size() == topology_.points().nb() );

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
AdaptationManager<type>::synchronize()
{
  // synchronize all the global point indices between processors
  // the root processor needs to know how many vertices there are in total
  if (rank_ == 0)
  {
    // receive the total number of vertices from all processors

    // receive the globals on the boundary interface

    // assign new global indices for the interface points

    // send global indices to all processors
  }
  else
  {
    // receive global indices for these points
    // and assign them
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
  for (index_t k=0;k<facets.nb();k++)
  {
    if (!facets.boundary(k)) continue;
    facets.retrieve(k,facet);
    if (BoundaryUtils::geometryFacet( topology_.points() , facet.data() , facet.size())!=nullptr) continue;

    for (index_t j=0;j<facet.size();j++)
      pts.push_back(facet[j]);
  }
  uniquify(pts);
  std::sort(pts.begin(),pts.end());

  // fix all the points
  for (index_t j=0;j<pts.size();j++)
    topology_.points().set_fixed(pts[j],true);

  // move all the fixed points to the beginning of the list
  // this should also adjust the global indices
  // but we'll need to adjust the metrics
  std::map<index_t,index_t> point_map;
  topology_.move_to_front( pts , &point_map );

  // map the metrics
  std::vector<VertexMetric> mapped_metrics(metrics_.size());
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

  bool done = false;
  for (index_t pass=0;pass<params_.max_passes();pass++)
  {
    // fix the boundary of the topology
    // and move the fixed points to the beginning of the points structure
    fix_boundary();

    // do the adaptation

    mpi::barrier();

    // synchronize all the global point indices


    mpi::barrier();

    // analyze whether we are done
    done = analyze();
    if (done) break;
  }

  // balance the mesh between the partitions
  balance();
}

template<typename type>
void
AdaptationManager<type>::retrieve( Topology<type>& topology ) const
{
  // fill the topology with the elements
  // fill the points with the vertices
  // assign the globals into the points
  avro_implement;
}


template class AdaptationManager<Simplex>;

#endif

} // avro

#endif
