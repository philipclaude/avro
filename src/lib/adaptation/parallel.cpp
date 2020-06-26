#include "common/mpi.hpp"
#include "common/process.h"

#include "adaptation/parallel.h"
#include "adaptation/parameters.h"

#include "element/simplex.h"

#include "library/meshb.h"

#include "mesh/facets.h"
#include "mesh/mpi_tags.h"
#include "mesh/partition.h"
#include "mesh/points.h"
#include "mesh/topology.h"

#include "numerics/matrix.h"

namespace avro
{

#define ACTIVE 1
#define INACTIVE 0
#define TAG_ACTIVE 100
#define TAG_REDUCE_PARTITION 200

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
  if (rank==0) printf("--> root: topology.nb() = %lu\n",topology_in.nb());
  else if (level > 0) avro_assert( topology_in.nb() == 0 ); // these are the empty interfaces

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

  repartition:

  // determine if an interface is needed
  std::vector<int> active( nb_rank , INACTIVE );
  if (topology_in.nb() <= partition_size )
  {
    // only turn on processor 1
    nb_partition = 1;
    if (topology_in.nb()>0)
      active[1] = ACTIVE;
  }
  else
  {
    #if 0
    index_t nb_partition0 = nb_partition;
    nb_partition = topology_in.nb()/partition_size;
    if (nb_partition>nb_partition0) nb_partition = nb_partition0;
    #endif

    // only turn on processors which are necessary
    for (index_t j=0;j<nb_partition;j++)
      active[j+1] = ACTIVE;
  }

  // communicate which processors are active
  if (rank==0)
  {
    for (index_t k=1;k<nb_rank;k++)
      mpi::send( mpi::blocking{} , active , k , TAG_ACTIVE );
  }
  else
    active = mpi::receive<std::vector<int>>( 0 , TAG_ACTIVE );
  mpi::barrier();

  // determine if we are done
  index_t nb_active = 0;
  for (index_t k=0;k<active.size();k++)
    nb_active += active[k];
  if (nb_active == 0) return 0; // we're done!
  printf("nb_active processors = %lu\n",nb_active);

  std::shared_ptr<AdaptationInterface<type> > interface = nullptr;
  std::shared_ptr<AdaptationChunk<type> > partition = nullptr;

  // initialize the interface
  // even non-root processors will maintain an interface so recursive
  // calls into this function can still be made
  interface = std::make_shared< AdaptationInterface<type> >(dim,udim,number);

  int reduce_partitions = 0;
  if (rank == 0)
  {
    avro_assert( topology_in.nb() > 0 );

    // copy in all the global points into the partion
    // we will prune unused points later
    topology_in.points().copy( interface->points() );

    std::vector< std::shared_ptr<Topology_Partition<type>> >parts( nb_partition );

    printf("--> partitioning level %lu with %lu elements into %d partitions\n",level,topology_in.nb(),nb_partition);
    Partition<type> partition(topology_in);
    partition.compute(nb_partition);

    // extract the partitions
    partition.get(parts);

    std::vector<std::set<index_t>> halo_points( nb_partition );
    partition.compute_interface_points( halo_points );

    for (index_t k=0;k<parts.size();k++)
    {
      parts[k]->set_entities(entities);

      // extract the halo for this partition and map to local indices
      std::vector<index_t> halo( halo_points[k].begin() , halo_points[k].end() );
      for (index_t j=0;j<halo.size();j++)
        halo[j] = parts[k]->global2local( halo[j] );

      // compute the set of elements touching a halo point
      parts[k]->neighbours().fromscratch() = true;
      parts[k]->neighbours().compute();
      parts[k]->inverse().build();
      std::vector<index_t> crust;
      parts[k]->compute_crust( halo , crust );

      // check if more than half of the partition is crust
      if (2*crust.size() > parts[k]->nb())
      {
        reduce_partitions = 1;
      }

      // add the crust to the interface
      // crust elements are in local indices, convert to global
      for (index_t j=0;j<crust.size();j++)
      {
        std::vector<index_t> s = parts[k]->get( crust[j] );
        for (index_t i=0;i<s.size();i++)
          s[i] = parts[k]->local2global(s[i]);
        interface->Topology<index_t>::add(s.data(),s.size());
      }
      printf("partition %lu has %lu elements with %lu elements in the crust\n",k+1,parts[k]->nb(),crust.size());

      // remove crust elements from the partition since they do not get adapted
      parts[k]->move_to_front( halo );
      parts[k]->remove_elements( crust );
      parts[k]->remove_unused();
    }

    // remove unused points in the interface (recall we used global indices when adding crust elements)
    interface->remove_unused();

    // send the partitions to each processor
    for (index_t k=0;k<parts.size();k++)
    {
      printf("--> sending partition %lu to processor %lu\n",k,k+1);
      parts[k]->send( comm , k+1 );
    }
  }
  else
  {
    // wait for the root to send us our partition
    partition = std::make_shared<AdaptationChunk<type>>(dim,udim,number);
    partition->set_entities(entities);
    if (active[rank]==ACTIVE)
    {
      partition->receive( comm , 0 );
      printf("--> received partition %lu with %lu elements\n",rank,partition->nb());
    }
    else
    {
      // nothing to receive since this processor is inactive
    }
  }
  mpi::barrier();

  if (rank==0)
  {
    // the root needs to inform the processors whether the partitions should be reduced
    for (index_t k=1;k<nb_rank;k++)
      mpi::send( mpi::blocking{} , reduce_partitions , k , TAG_REDUCE_PARTITION );
  }
  else
  {
    reduce_partitions = mpi::receive<int>( 0 , TAG_REDUCE_PARTITION );
  }
  mpi::barrier();

  if (reduce_partitions==1)
  {
    printf("*** reducing number of partitions ***\n");
    nb_partition--;
    goto repartition;
  }
  mpi::barrier();


  // do the serial adaptation
  if (rank == 0)
  {
    // receive the adapted partitions from the processors
    // note: the topologies are all in local indexing
    std::vector< std::shared_ptr<Topology_Partition<type>> > parts( nb_partition );
    for (index_t k=0;k<nb_partition;k++)
    {
      parts[k] = std::make_shared<Topology_Partition<type>>(dim,udim,number);
      parts[k]->set_entities(entities);
      parts[k]->receive( comm , k+1 );

      printf("--> received adapted partition from processor %lu\n",k+1);

      // compute the mantle and add it to the interface
      #if 0
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
      #endif

      // add the mantle elements to the interface
      // TODO
    }
  }
  else
  {
    // determine the fixed points and move them to the front
    // be sure to adjust local2global maps
    if (active[rank]==ACTIVE)
    {
      // adapt the mesh
      partition->adapt();

      // send the mesh to the root processor
      printf("--> sending adapted partition %lu with %lu elements back to root\n",rank,partition->nb());
      partition->send( comm , 0 );

      // extract all elements on the interface and then all vertices in this set of elements

      // send the set of interface elements as those in the ball of any vertex we found earlier
      //mpi::send(); // TODO
    }
    else
    {
      avro_assert( partition->nb() == 0 );
    }
  }
  mpi::barrier();

  // adapt the interface in parallel
  Points interface_points;
  Topology<type> interface_out( interface_points , interface->number() );

  interface->neighbours().fromscratch() = true;
  interface->neighbours().compute();

  // wait for all processors to make the recursive call
  mpi::barrier();

  if (rank==0)
  {
    library::meshb out;
    Mesh mesh(number,dim);
    interface->points().copy(mesh.points());
    mesh.add(interface);
    //out.write(mesh,"interface-"+std::to_string(level)+".mesh",false);
  }

  adaptp( *interface , metrics , params , interface_out , level+1 );

  // accumulate the adapted interface into the global output opology
  // TODO

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
