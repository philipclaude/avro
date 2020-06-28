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

#define ACTIVE 1
#define INACTIVE 0

#ifdef AVRO_MPI

template<typename type>
class AdaptationChunk : public Topology_Partition<type>
{
public:
  AdaptationChunk( coord_t dim , coord_t udim , coord_t number ) :
    Topology_Partition<type>(dim,udim,number)
  {}

  void fix_interface()
  {
    Facets facets(*this);
    facets.compute();

    //std::vector<index_t> pts;
    std::vector<index_t> facet(this->number());
    for (index_t k=0;k<facets.nb();k++)
    {
      if (!facets.boundary(k)) continue;
      facets.retrieve(k,facet);
      Entity* entity = BoundaryUtils::geometryFacet(this->points(),facet.data(),facet.size());
      if (entity!=nullptr) continue;
      if (fixed_facet(facet,this->points())) continue;
      for (index_t j=0;j<facet.size();j++)
        fixed_points_.push_back(facet[j]);
    }
    uniquify(fixed_points_);

    // move the boundary points to the front..they will be fixed during the adaptation
    // they need to be in the front so that the local2global indices are not touched
    // which would occur if they were placed after a vertex that gets collapsed
    std::map<index_t,index_t> idx;
    this->move_to_front(fixed_points_,&idx);
    avro_assert( this->all_points_accounted() );

    // recompute the facets...is there a way to not do this twice?
    facets.compute();
    fixed_points_.clear();
    for (index_t k=0;k<facets.nb();k++)
    {
      if (!facets.boundary(k)) continue;
      facets.retrieve(k,facet);
      Entity* entity = BoundaryUtils::geometryFacet(this->points(),facet.data(),facet.size());
      if (entity!=nullptr) continue;
      if (fixed_facet(facet,this->points())) continue;
      for (index_t j=0;j<facet.size();j++)
        fixed_points_.push_back(facet[j]);
    }
    uniquify(fixed_points_);

    // all of these boundary points should already be at the beginning of the point list
    for (index_t k=0;k<fixed_points_.size();k++)
    {
      avro_assert( fixed_points_[k] == k ); // ensures boundary points are at the beginning
      this->points().set_fixed(k,true);
    }

    // adjust the local2global indices
    std::vector<index_t> global = this->local2global_;
    for (index_t k=0;k<this->points().nb();k++)
    {
      this->local2global_[idx.at(k)] = global[k];
      this->global2local_[global[k]] = idx.at(k);
    }
  }

  void adapt()
  {
    // call the serial adaptation

  }

  void unfix_interface()
  {
    for (index_t k=0;k<this->points().nb();k++)
      this->points().set_fixed(k,true);
    for (index_t k=0;k<fixed_points_.size();k++)
      this->points().set_fixed(fixed_points_[k],false);
  }

  Entity* lookup( index_t identifier ) const
  {
    avro_implement;
    return nullptr;
  }

  std::vector<VertexMetric>& metric() { return metrics_; }

private:
  std::vector<VertexMetric> metrics_;
  std::vector<index_t> fixed_points_;
};

int
estimate_partition_size( int nb_partition )
{
  return 1000; // TODO analyze disk space or requested memory from user to determine how many elements
}

template<typename type>
int
adaptp( Topology<type>& topology_in , const std::vector<VertexMetric>& metrics , AdaptationParameters& params , Topology<type>& topology_out , index_t level , int finished0 )
{
  index_t rank = mpi::rank();
  mpi::communicator& comm = ProcessMPI::get_comm();
  index_t nb_rank = mpi::size();
  if (rank==0) printf("--> root: topology.nb() = %lu, finished = %d\n",topology_in.nb(),finished0);
  else if (level > 0) avro_assert( topology_in.nb() == 0 ); // these are the empty interfaces

  const coord_t dim = topology_in.points().dim();
  const coord_t udim = topology_in.points().udim();
  const coord_t number = topology_in.number();

  //if (rank==0) topology_in.points().print(true);

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

  // sigh...I know goto's are not good programming practice
  // but it's just so convenient here
  // this is the goto used in case the number of partitions needs to be reduced
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
      mpi::send( mpi::blocking{} , active , k , TAG_MISC );
  }
  else
    active = mpi::receive<std::vector<int>>( 0 , TAG_MISC );
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

  std::vector<std::set<index_t>> halo_points( nb_partition );
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

    // extract the partitions along with the interface points
    partition.get(parts);
    partition.compute_interface_points( halo_points );

    for (index_t k=0;k<parts.size();k++)
    {
      parts[k]->set_entities(entities);

      // extract the halo for this partition and map to local indices
      #if 0
      std::vector<index_t> halo( halo_points[k].begin() , halo_points[k].end() );
      for (index_t j=0;j<halo.size();j++)
        halo[j] = parts[k]->global2local( halo[j] );
      #else
      std::vector<index_t> halo;
      std::set<index_t>::iterator it;
      for (it=halo_points[k].begin();it!=halo_points[k].end();it++)
      {
        if (parts[k]->points().fixed( *it ) ) continue;
        halo.push_back( parts[k]->global2local( *it ) );
      }
      #endif

      // compute the set of elements touching a halo point
      parts[k]->neighbours().fromscratch() = true;
      parts[k]->neighbours().compute();
      parts[k]->inverse().build();
      std::vector<index_t> crust;
      parts[k]->compute_crust( halo , crust );

      // check if more than half of the partition is crust
      //if (2*crust.size() > parts[k]->nb())
      if (crust.size() >= parts[k]->nb())
      {
        // initiate the flag that partitions need to be reduced
        reduce_partitions = 1;
      }

      // add the crust to the interface
      // crust elements are in local indices, convert to global
      for (index_t j=0;j<crust.size();j++)
      {
        std::vector<index_t> s = parts[k]->get( crust[j] );
        for (index_t i=0;i<s.size();i++)
          s[i] = parts[k]->local2global(s[i]);
        //if (level==0)
        interface->Topology<index_t>::add(s.data(),s.size());
      }
      printf("partition %lu has %lu elements with %lu elements in the crust\n",k+1,parts[k]->nb(),crust.size());

      // remove crust elements from the partition since they do not get adapted
      std::map<index_t,index_t> idx;
      parts[k]->move_to_front( halo , &idx );
      parts[k]->map_indices(idx);
      parts[k]->remove_elements( crust );
      std::vector<index_t> pts;
      parts[k]->remove_unused(&pts);
      parts[k]->remove_indices(pts);

      bool check = parts[k]->check( interface->points() );
      avro_assert( check );
    }

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
      mpi::send( mpi::blocking{} , reduce_partitions , k , TAG_MISC );
  }
  else
  {
    // wait for the root to tell us if we need to go back
    reduce_partitions = mpi::receive<int>( 0 , TAG_MISC );
  }
  mpi::barrier();

  // check if partitions need to be reduced and go back to the top!
  if (reduce_partitions==1)
  {
    printf("*** reducing number of partitions ***\n");
    nb_partition--;
    goto repartition;
  }
  mpi::barrier();

  // do the serial adaptation
  int finished = finished0;
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

      bool check = parts[k]->check( interface->points() );
      avro_assert( check );

      printf("--> received adapted partition with %lu elements from processor %lu\n",parts[k]->nb(),k+1);

      if (finished0==0)
      {
        // all points were unfixed before they were sent here


        // compute the mantle and add it to the interface
        // the mantle consists of all elements touching a vertex on the adapted partition boundary
        // in which the boundary is identified from any facet that is on the actual geometry boundary
        // and which is not currently fixed
        parts[k]->neighbours().fromscratch() = true;
        parts[k]->neighbours().compute();
        parts[k]->inverse().build();
        std::vector<index_t> interior, exterior, crust;
        parts[k]->compute_mantle( interior , exterior , crust , halo_points[k] );

        //parts[k]->points().print(true);

        printf("--> there are %lu elements in the crust with %lu interior points and %lu exterior points\n",crust.size(),interior.size(),exterior.size());
        if (2*crust.size() > parts[k]->nb() && finished0==0)
        {
          printf("*** crust is too big --> will perform next adaptation in serial ***\n");
          finished = 1;
        }

        std::map<index_t,index_t> local2global;
        for (index_t j=0;j<exterior.size();j++)
          local2global.insert( {exterior[j],parts[k]->local2global(exterior[j])} );

        // add all the points in the interior to the interior mesh
        // exterior points should already be added
        for (index_t j=0;j<interior.size();j++)
        {
          index_t q = interface->points().nb();
          index_t p = interior[j];
          interface->points().create( parts[k]->points()[p] );
          interface->points().set_entity( q , parts[k]->points().entity(p) );
          interface->points().set_param( q , parts[k]->points().u(p) );
          interface->points().set_fixed( q , parts[k]->points().fixed(p) );
          local2global.insert( {p,q} );
        }

        for (index_t j=0;j<crust.size();j++)
        {
          // retrieve the mantle element and add it to the interface
          std::vector<index_t> s = parts[k]->get(crust[j]);

          // the mantle element s is in local partition indexing...
          // it needs to be converted to interface coordinates
          // which are "global" (until we remove unused points below)
          for (index_t i=0;i<s.size();i++)
            s[i] = local2global.at( s[i] );

          // add the mantle element to the interface
          interface->add( s.data() , s.size() );
        }

        // remove the mantle from the adapted partition (these elements will be added by the interface)
        parts[k]->remove_elements( crust );

        // add the shrunk partition to the global topology
        printf("[TODO]: add shrunk partition to the global topology\n");
      }
      else
      {
        avro_assert( nb_partition==1 );

        // we were not requested an interface mesh

        // add this adapted partition to the global topology
        printf("[TODO]: add entire single adapted partition to output topology!\n");
      }

    }
  }
  else
  {
    // determine the fixed points and move them to the front
    // be sure to adjust local2global maps
    if (active[rank]==ACTIVE)
    {
      // the fixed points should have already been moved to the front
      partition->fix_interface();

      // now we can retrieve the metrics from the global list
      // TODO

      // adapt the mesh
      partition->adapt();

      // unfix the partition
      for (index_t j=0;j<partition->points().nb();j++)
        partition->points().set_fixed(j,false);

      partition->unfix_interface();

      // send the mesh to the root processor
      printf("--> sending adapted partition %lu with %lu elements back to root\n",rank,partition->nb());
      partition->send( comm , 0 );
    }
    else
    {
      avro_assert( partition->nb() == 0 );
    }
  }
  mpi::barrier();

  if (rank==0)
  {
    // the root needs to inform the processors whether the partitions should be reduced
    for (index_t k=1;k<nb_rank;k++)
      mpi::send( mpi::blocking{} , finished , k , TAG_MISC );
  }
  else
  {
    // wait for the root to tell us if we need to go back
    finished = mpi::receive<int>( 0 , TAG_MISC );
  }
  mpi::barrier();

  // adapt the interface in parallel
  Points interface_points;
  Topology<type> interface_out( interface_points , interface->number() );

  // remove unused points in the interface (recall we used global indices when adding crust elements)
  interface->remove_unused();
  interface->neighbours().fromscratch() = true;
  interface->neighbours().compute();

  if (rank==0)
  {
    library::meshb out;
    Mesh mesh(number,dim);
    interface->points().copy(mesh.points());
    mesh.add(interface);
    out.write(mesh,"interface-"+std::to_string(level)+".mesh",false);
  }

  // retrieve the metrics passed in to the interface adaptation
  // TODO

  for (index_t k=0;k<interface->points().nb();k++)
    interface->points().set_fixed(k,false);

  // assign the boundary of the interface as fixed so that it will not be considered a partition "boundary"
  Facets facets(*interface);
  facets.compute();
  std::vector<index_t> facet(number);
  index_t nb_bnd = 0;
  std::vector<index_t> bnd_pts;
  for (index_t k=0;k<facets.nb();k++)
  {
    if (!facets.boundary(k)) continue;
    facets.retrieve(k,facet);
    nb_bnd++;
    for (index_t j=0;j<facet.size();j++)
    {
      bnd_pts.push_back(facet[j]);
      interface->points().set_fixed( facet[j] , true );
    }
  }
  uniquify(bnd_pts);
  if (rank==0)
  {
    printf("number of boundary facets = %lu out of %lu\n",nb_bnd,facets.nb());
    printf("number of boundary points = %lu out of %lu\n",bnd_pts.size(),interface->points().nb());
  }

  // if the finished flag was initiated then the next adaptation must occur in serial
  if (finished==1)
    params.nb_partition() = 1;

  // we started the adaptation assuming the interface needed to be adapted
  if (finished0==0)
  {
    // adapt the interface
    adaptp( *interface , metrics , params , interface_out , level+1 , finished );

    /*if (rank==0)
      interface->points().print(true);
    else
      avro_assert( interface->points().nb()==0 );*/

    // accumulate the adapted interface into the global output topology
    printf("[TODO]: add adapted interface to global topology\n");
  }
  else
  {
    // the interface should have been completely added by parts[0] (above)
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
template int adaptp( Topology<Simplex>& , const std::vector<VertexMetric>& , AdaptationParameters& , Topology<Simplex>& , index_t level , int );

#endif

} // avro
