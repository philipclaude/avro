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
class AdaptationChunk : public Topology_Partition<type>
{
public:
  AdaptationChunk( coord_t dim , coord_t udim , coord_t number ) :
    Topology_Partition<type>(dim,udim,number),
    points_(dim,udim)
  {}

  void adapt()
  {
    // call the serial adaptation

  }

  Entity* lookup( index_t identifier ) const
  {
    avro_implement;
    return nullptr;
  }

  std::vector<VertexMetric>& metric() { return metrics_; }

private:
  Points points_;
  std::vector<index_t> local2global_;
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
  printf("rank = %lu, nb_rank = %lu\n",rank,nb_rank);

  const coord_t dim = topology_in.points().dim();
  const coord_t udim = topology_in.points().udim();
  const coord_t number = topology_in.number();

  // determine the geometric entities
  std::set<Entity*> entities0;
  for (index_t k=0;k<topology_in.points().nb();k++)
  {
    if (topology_in.points().entity(k)==nullptr) continue;
    entities0.insert( topology_in.points().entity(k) );
  }
  std::vector<Entity*> entities(entities0.begin(),entities0.end());

  // determine if we need to partition the mesh
  int nb_partition = params.nb_partition();
  int partition_size = params.partition_size();
  if (partition_size<0)
    partition_size = estimate_partition_size(nb_partition);
  printf("partition size = %d\n",partition_size);

  // determine if an interface is needed
  bool interface_needed = true;
  if (topology_in.nb() <= partition_size )
  {
    interface_needed = false;
    nb_partition     = 1;
  }

  std::shared_ptr<AdaptationInterface<type> > interface = nullptr;
  std::shared_ptr<AdaptationChunk<type> > partition = nullptr;
  if (rank == 0)
  {
    // partition the mesh
    printf("partitioning level %lu with %lu elements into %d partitions\n",level,topology_in.nb(),nb_partition);
    Partition<type> partition(topology_in);
    partition.compute(nb_partition);
    printf("done\n");

    // initialize the interface
    interface = std::make_shared< AdaptationInterface<type> >(dim,udim,number);
    topology_in.points().copy( interface->points() );

    // extract the partitions
    std::vector< std::shared_ptr<Topology_Partition<type>> >parts( nb_partition );
    partition.get(parts);

    for (index_t k=0;k<parts.size();k++)
    {
      parts[k]->set_entities(entities);

      // compute the crust of each partition
      // TODO, need a better method for determining partition boundaries without depending on the geometry
      // implement this as a method of the global partition class
      parts[k]->compute_crust();

      const std::vector<index_t>& crust = parts[k]->crust();

      // add the crust to the interface
      // everything is still in global coordinates
      index_t np = interface->points().nb();
      for (index_t j=0;j<crust.size();j++)
      {
        std::vector<index_t> s = parts[k]->get( crust[j] );
        /*for (index_t i=0;i<s.size();i++)
          s[i] = interface->add_point( s[i] , parts[k]->points() , np );*/

        for (index_t i=0;i<s.size();i++)
          s[i] = parts[k]->local2global(s[i]);
        interface->Topology<index_t>::add(s.data(),s.size());
      }

      // remove crust elements from the partition since they do not get adapted

      print_inline( parts[k]->halo() , "halo" );
      parts[k]->move_to_front( parts[k]->halo() );
      parts[k]->remove_elements( crust );
      parts[k]->remove_unused();
    }

    // send the partitions to each processor
    printf("sending partitions to processors...\n");
    for (index_t k=1;k<nb_rank;k++)
      parts[k-1]->send( comm , k );
    printf("done\n");
  }
  else
  {
    // wait for the root to send us our partition
    printf("waiting to receive partitions on rank %lu\n",rank);
    partition = std::make_shared<AdaptationChunk<type>>(dim,udim,number);
    partition->set_entities(entities);
    partition->receive( comm , 0 );
  }
  mpi::barrier();

  // do the serial adaptation
  if (rank == 0)
  {
    // receive the adapted partitions from the processors
    // note: the topologies are all in local indexing
    std::vector< std::shared_ptr<Topology_Partition<type>> > parts( nb_partition );
    for (index_t k=0;k<nb_rank-1;k++)
    {
      parts[k] = std::make_shared<Topology_Partition<type>>(dim,udim,number);
      parts[k]->set_entities(entities);
      parts[k]->receive( comm , k+1 );

      printf("received adapted partition from processor %lu\n",k+1);

      // compute the mantle and add it to the interface
      parts[k]->compute_mantle();

      const std::vector<index_t>& mantle = parts[k]->mantle();

      for (index_t j=0;j<mantle.size();j++)
      {
        // retrieve the mantle element and add it to the interface
        std::vector<index_t> s = parts[k]->get(mantle[j]);

        // the mantle element s is in local partition coordinates
        // it needs to be converted to local interface coordinates
        for (index_t i=0;i<s.size();i++)
        {
          // determine if s[i] is a
        }

      }

      // remove the mantle from the adapted partition (it will be added by the interface)
      parts[k]->remove_elements( mantle );

      // add the mantle elements to the interface
      // TODO
    }


  }
  else
  {
    // determine the fixed points and move them to the front
    // be sure to adjust local2global maps

    // adapt the mesh
    partition->adapt();

    // send the mesh to the root processor
    printf("sending adapted partition %lu back to root\n",rank);
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

    printf("adapting interface:\n");
    interface->remove_unused();
    interface->Table<index_t>::print();
    interface->neighbours().fromscratch() = true;
    interface->neighbours().compute();

    adaptp( *interface , metrics , params , interface_out , level+1 );

    // accumulate the adapted interface into the global output opology
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
