#include "common/mpi.hpp"
#include "common/process.h"

#include "adaptation/parallel.h"
#include "adaptation/parameters.h"

#include "element/simplex.h"

#include "mesh/mpi_tags.h"
#include "mesh/partition.h"
#include "mesh/points.h"
#include "mesh/topology.h"

#include "numerics/matrix.h"

namespace avro
{

#ifdef AVRO_MPI

template<typename type>
class PartitionInterface : public Topology<type>
{
public:
  PartitionInterface( const Topology<type>& topology ) :
    Topology<type>(points_,topology.number()),
    points_(topology.points().dim()),
    topology_(topology)
  {}

  void add( );

private:
  Points points_;
  const Topology<type>& topology_;
  std::map<index_t,index_t> local2global_;
};

template<typename type>
class PartitionChunk : public Topology_Partition<type>
{
public:
  PartitionChunk( const Topology<type>& topology ) :
    Topology_Partition<type>(points_.dim(),points_.udim(),topology.number()),
    points_(topology.points().dim()),
    topology_(topology),
    boundary_(points_,topology.number()-1)
  {
    initialize();
  }

  void initialize()
  {
    // add all the elements from the topology

    // create the points and the local2global indices for the points on the boundary

    // extract the boundary
  }

  void fix_points( const PartitionInterface<type>& interface )
  {
    // fix all points in the interface that are points in this topology
  }

  void extract_metric( std::vector<VertexMetric>& global_metrics )
  {
    // pick out the vertex metrics for this topology from the global list
  }

  void adapt()
  {
    // call the serial adaptation

  }

  Entity* lookup( index_t identifier ) const
  {
    avro_implement;
    return nullptr;
  }

  const Topology<type>& boundary() const { return boundary_; }
  std::vector<VertexMetric>& metric() { return metrics_; }

private:
  Points points_;
  const Topology<type>& topology_;
  std::vector<index_t> local2global_;
  Topology<type> boundary_;
  std::vector<VertexMetric> metrics_;
};

int
estimate_partition_size( int nb_partition )
{
  return 1000; // TODO analyze disk space or requested memory from user to determine how many elements
}

template<typename type>
int
adaptp( Topology<type>& topology_in , const std::vector<VertexMetric>& metrics , AdaptationParameters& params , Topology<type>& topology_out , index_t level )
{
  index_t rank = mpi::rank();
  mpi::communicator& comm = ProcessMPI::get_comm();
  index_t nb_rank = mpi::size();
  printf("rank = %lu\n",rank);

  // determine if we need to partition the mesh
  int nb_partition = params.nb_partition();
  int partition_size = params.partition_size();
  if (partition_size<0)
    partition_size = estimate_partition_size(nb_partition);

  // determine if an interface is needed
  bool interface_needed = true;
  if (topology_in.nb() <= partition_size )
  {
    interface_needed = false;
    nb_partition     = 1;
  }

  std::shared_ptr<PartitionChunk<type>> partition = nullptr;
  if (rank == 0)
  {
    // partition the mesh
    printf("partitioning level %lu with %lu partitions\n",level,nb_partition);
    Partition<type> partition(topology_in);
    partition.compute(nb_partition);
    printf("done\n");

    // extract the partitions
    std::vector< std::shared_ptr<Topology_Partition<type>> >parts( nb_partition );
    partition.get(parts);

    // send the partitions to each processor
    printf("sending partitions to processors...\n");
    for (index_t k=1;k<nb_rank;k++)
      parts[k]->send( comm , k );
    printf("done\n");
  }
  else
  {
    // wait for the root to send us our partition
    printf("waiting to receive partitions on rank %lu\n",rank);
    partition = std::make_shared<PartitionChunk<type>>(topology_in);
    partition->receive( comm , 0 );
  }
  mpi::barrier();
  return 1;

  // do the serial adaptation
  std::shared_ptr<PartitionInterface<type> > interface = nullptr;
  if (rank == 0)
  {
    // wait for the result
    for (index_t k=1;k<nb_rank;k++)
    {
      Topology<type> tk( topology_out.points() , topology_in.number() );
      tk.receive( comm , k );

      // append the result to the outgoing topology
      for (index_t j=0;j<tk.nb();j++)
      {
        // TODO (remember to offset the indices by the current number of points)
      }

      // append the vertices too
      for (index_t j=0;j<tk.points().nb();j++)
      {
        // TODO
      }
    }

    if (interface_needed)
    {
      // initialize the interface
      interface = std::make_shared< PartitionInterface<type> >(topology_out);

      // receive the set of interface elements from each processor
      //std::vector<index_t> interface_elems = mpi::recv();

      // add the interface elements to the interface
      // TODO

    }
  }
  else
  {
    // adapt the mesh
    partition->adapt();

    // send the mesh to the root processor
    partition->send( comm , 0 );

    // extract all elements on the interface and then all vertices in this set of elements

    // send the set of interface elements as those in the ball of any vertex we found earlier
    //mpi::send(); // TODO
  }
  mpi::barrier();

  if (rank == 0 && interface_needed)
  {
    // adapt the interface in parallel
    Points interface_points;
    Topology<type> interface_out( interface_points , interface->number() );
    adaptp( *interface , metrics , params , interface_out , level+1 );

    // accumulate the interface into the global topology
    // TODO
  }
  mpi::barrier();

  return 0;
}

template<typename type>
AdaptationManager<type>::AdaptationManager( Topology<type>& topology ,
  std::vector<VertexMetric>& metrics ,
  AdaptationParameters& params) :
  topology_(topology),
  metrics_(metrics),
  params_(params)
{

}


template class AdaptationManager<Simplex>;
template int adaptp( Topology<Simplex>& , const std::vector<VertexMetric>& , AdaptationParameters& , Topology<Simplex>& , index_t level );

#endif

} // avro
