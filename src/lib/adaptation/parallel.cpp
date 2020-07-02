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
  if (rank == 0)
    printf("*** nb_active processors = %lu ***\n",nb_active);

  std::shared_ptr<AdaptationInterface<type> > interface = nullptr;
  std::shared_ptr<AdaptationChunk<type> > partition = nullptr;

  // initialize the interface
  // even non-root processors will maintain an interface so recursive
  // calls into this function can still be made
  interface = std::make_shared< AdaptationInterface<type> >(dim,udim,number);

  std::vector<std::set<index_t>> partition_boundary_points( nb_partition );
  std::set<index_t> boundary_points;
  int reduce_partitions = 0;
  if (rank == 0)
  {
    avro_assert( topology_in.nb() > 0 );

    // copy in all the global points into the partion
    // we will prune unused points later
    //topology_in.points().copy( interface->points() );

    printf("--> partitioning level %lu with %lu elements into %d partitions\n",level,topology_in.nb(),nb_partition);
    std::vector< std::shared_ptr<Topology_Partition<type>> >parts( nb_partition );
    Partition<type> partition(topology_in);
    partition.compute(nb_partition);

    // extract the partitions along with the points on the partition boundaries
    partition.get(parts);
    partition.compute_interface_points( partition_boundary_points );

    for (index_t k=0;k<parts.size();k++)
    {
      parts[k]->set_entities(entities);

      // extract the halo for this partition and map to local indices
      #if 0
      std::vector<index_t> halo( partition_boundary_points[k].begin() , partition_boundary_points[k].end() );
      for (index_t j=0;j<halo.size();j++)
      {
        boundary_points.insert( halo[j] );
        halo[j] = parts[k]->global2local( halo[j] );
      }
      #else
      std::vector<index_t> halo;
      std::set<index_t>::iterator it;
      for (it=partition_boundary_points[k].begin();it!=partition_boundary_points[k].end();it++)
      {
        if (parts[k]->points().fixed( *it ) ) continue;
        halo.push_back( parts[k]->global2local( *it ) );
        boundary_points.insert( *it );
      }
      #endif

      // move the partition boundary points to the front
      std::map<index_t,index_t> idx;
      std::vector<index_t> pts;
      parts[k]->move_to_front( halo , &idx );
      parts[k]->map_indices(idx); // maps the local2global indices based on the shift to the front
      parts[k]->remove_unused(&pts);  // removes extra points
      parts[k]->remove_indices(pts);  // removes local2global indices associated with extra points

      bool check = parts[k]->check( topology_in.points() );
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

  // communicate whether partitions should be reduces
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

  std::map<index_t,index_t> output2interface;
  std::map<index_t,index_t> interface2output;

  // do the serial adaptation
  int finished = finished0;
  if (rank == 0)
  {

    // add the interface points to the global topology
    std::map<index_t,index_t> global2output;
    for (std::set<index_t>::iterator it=boundary_points.begin();it!=boundary_points.end();++it)
    {
      index_t p = *it;
      index_t q = topology_out.points().nb();
      topology_out.points().create( topology_in.points()[p] );
      topology_out.points().set_entity( q , topology_in.points().entity(p) );
      topology_out.points().set_param( q , topology_in.points().u(p) );
      topology_out.points().set_fixed( q , topology_in.points().fixed(p) );
      global2output.insert( {p,q} );
    }
    printf("added %lu points to output topology\n",global2output.size());

    std::vector< std::shared_ptr<Topology_Partition<type>> > parts( nb_partition );
    for (index_t k=0;k<nb_partition;k++)
    {
      // receive the adapted partition from the processor
      // note: the topology is in local indexing
      // with fixed points in the beginning of the points container
      parts[k] = std::make_shared<Topology_Partition<type>>(dim,udim,number);
      parts[k]->set_entities(entities);
      parts[k]->receive( comm , k+1 );

      // check the local2global indexing is correct
      bool check = parts[k]->check( topology_in.points() );
      avro_assert( check );

      printf("--> received adapted partition with %lu elements from processor %lu\n",parts[k]->nb(),k+1);

      // add the vertices to the output topology and
      // compute the map between the local partition indices to the output indices
      std::map<index_t,index_t> part2output;
      for (index_t j=0;j<parts[k]->points().nb();j++)
      {
        index_t p = parts[k]->local2global(j);

        // this vertex already exists
        if (global2output.find(p)!=global2output.end())
        {
          part2output.insert( {j,global2output.at(p)} );
          continue;
        }

        // vertex does not yet exist so add it
        index_t q = topology_out.points().nb();
        topology_out.points().create( parts[k]->points()[j] );
        topology_out.points().set_entity( q , parts[k]->points().entity(j) );
        topology_out.points().set_param( q , parts[k]->points().u(j) );
        topology_out.points().set_fixed( q , parts[k]->points().fixed(j) );
        part2output.insert( {j,q} );
      }

      // add all elements in the partition to the output topology
      // we will remove interface elements later
      for (index_t j=0;j<parts[k]->nb();j++)
      {
        // retrieve the element in the partition
        std::vector<index_t> s = parts[k]->get(j);

        // map the element indices to the output indices
        for (index_t i=0;i<s.size();i++)
          s[i] = part2output.at(s[i]);

        // add the element to the output topology
        topology_out.add( s.data() , s.size() );
      }
    } // loop over partitions

    // build the data structures for the output topology
    topology_out.build_structures();

    // identify all elements in the output topology that are:
    //   1) in the ball of the interface points
    //   2) in the ball of any point in an element in the ball of (1)
    std::set<index_t> crust;
    std::vector<index_t> ball1,ball2;
    for (std::set<index_t>::iterator it=boundary_points.begin();it!=boundary_points.end();++it)
    {
      index_t p = global2output.at( *it );
      ball1.clear();
      topology_out.inverse().ball( p , ball1 );
      for (index_t i=0;i<ball1.size();i++)
      {
        for (index_t j=0;j<topology_out.nv(ball1[i]);j++)
        {
          index_t q = topology_out(ball1[i],j);
          ball2.clear();
          topology_out.inverse().ball( q , ball2 );
          for (index_t k=0;k<ball2.size();k++)
            crust.insert(ball2[k]);
        }
      }
    }
    printf("there are %lu elements in the crust\n",crust.size());

    // add every element in the crust to the interface
    for (std::set<index_t>::iterator it=crust.begin();it!=crust.end();++it)
    {
      std::vector<index_t> s = topology_out.get( *it );

      // remap the indices to local interface indices
      for (index_t i=0;i<s.size();i++)
      {
        index_t p = s[i];
        if (output2interface.find(p) == output2interface.end())
        {
          // add the point to the interface
          index_t q = interface->points().nb();
          interface->points().create( topology_out.points()[p] );
          interface->points().set_entity( q , topology_out.points().entity(p) );
          interface->points().set_param( q , topology_out.points().u(p) );
          interface->points().set_fixed( q , topology_out.points().fixed(p) );
          output2interface.insert( {p,q} );
          interface2output.insert( {q,p} );
        }
        s[i] = output2interface.at(p);
      }
      interface->add( s.data() , s.size() );
    }
    avro_assert( interface->points().nb() == output2interface.size() );

    // remove the crust elements
    std::vector<index_t> crust0(crust.begin(),crust.end());
    topology_out.remove_elements(crust0);

    std::vector<index_t> pts; // the indices of the removed points (in the output)
    topology_out.determine_unused(pts);

    // check the original map
    for (std::map<index_t,index_t>::iterator it=interface2output.begin();it!=interface2output.end();++it)
    {
      index_t p_i = it->first;
      index_t p_o = it->second;
      real_t d0 = numerics::distance( interface->points()[p_i] , topology_out.points()[ p_o ] , dim );
      avro_assert( d0 < 1e-3 );
    }

    // copy the points so we can check the map
    Points P(dim);
    topology_out.points().copy(P);

    std::map<index_t,index_t> output2output;
    topology_out.move_to_front(pts,&output2output);

    // make sure the map is correct!
    for (std::map<index_t,index_t>::const_iterator it=output2output.begin();it!=output2output.end();++it)
    {
      index_t p0 = it->first;
      index_t p1 = it->second;

      real_t d = numerics::distance( P[p0] , topology_out.points()[p1] , dim );
      avro_assert( d<1e-3 );
    }
    topology_out.remove_points(pts);

    std::map<index_t,index_t> interface2output2;
    for (std::map<index_t,index_t>::const_iterator it=interface2output.begin();it!=interface2output.end();it++)
    {
      index_t p_i = it->first;
      index_t p_o = it->second;

      // compute the new output point
      index_t p_o2 = output2output.at(p_o);
      if (p_o2<pts.size()) continue;
      p_o2 -= pts.size();

      real_t d = numerics::distance( interface->points()[p_i] , topology_out.points()[p_o2] , dim );
      if (d > 1e-3)
      {
        interface->points().print(p_i);
        topology_out.points().print(p_o2);
      }
      interface2output2.insert({p_i,p_o2});
    }
    avro_assert( interface2output.size() == interface2output2.size() + pts.size() );
    interface2output = interface2output2;
    avro_assert( topology_out.all_points_accounted() );
  }
  else
  {
    // determine the fixed points and move them to the front
    // be sure to adjust local2global maps
    if (active[rank]==ACTIVE)
    {
      // the fixed points should have already been moved to the front
      //partition->fix_interface();

      // now we can retrieve the metrics from the global list
      // TODO

      // adapt the mesh
      partition->adapt();

      // unfix the partition before sending it to the root
      // everything in the partition will be fixed except for the partition
      // boundary, which is needed to extract the correct set of mantle elements
      //partition->unfix_interface();

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
  Points interface_points(dim,udim);
  Topology<type> interface_out( interface_points , interface->number() );

  // assign the boundary of the interface as fixed so that it will not be considered a partition "boundary"
  interface->build_structures();
  #if 0
  Facets facets(*interface);
  facets.compute();
  std::vector<index_t> facet(number);
  for (index_t k=0;k<facets.nb();k++)
  {
    if (!facets.boundary(k)) continue;
    facets.retrieve(k,facet);
    if (BoundaryUtils::geometryFacet( interface->points() , facet.data() , facet.size() )!=nullptr) continue;
    for (index_t j=0;j<facet.size();j++)
      interface->points().set_fixed( facet[j] , true );
  }

  std::vector<index_t> interface_fixed_points;
  for (index_t j=0;j<interface->points().nb();j++)
  {
    if (interface->points().fixed(j))
      interface_fixed_points.push_back(j);
  }
  interface->move_to_front( interface_fixed_points );
  #endif
  avro_assert( interface->all_points_accounted() );

  if (rank==0)
  {
    library::meshb out;
    Mesh mesh(number,dim);
    interface->points().copy(mesh.points());
    mesh.add(interface);
    out.write(mesh,"interface-"+std::to_string(level)+".mesh",false);

    topology_out.add_child( interface );
  }

  // make sure the interface2output map is correct before adapting the interface
  index_t nerror = 0;
  index_t ncheck = 0;
  for (index_t k=0;k<interface->points().nb();k++)
  {
    if (interface2output.find(k)==interface2output.end()) continue;
    ncheck++;

    index_t p = interface2output.at(k);
    real_t d = numerics::distance( interface->points()[k] , topology_out.points()[p] , dim );
    if (d>1e-3)
    {
      interface->points().print(k);
      topology_out.points().print(p);
      printf("distance = %g\n",d);
      nerror++;
    }
  }
  if (rank == 0 )
    printf("there were %lu errors out of %lu\n",nerror,ncheck);
  avro_assert( nerror == 0 );

  // retrieve the metrics passed in to the interface adaptation
  // TODO

  // if the finished flag was initiated then the next adaptation must occur in serial
  if (finished==1)
    params.nb_partition() = 1;

  // we started the adaptation assuming the interface needed to be adapted
  if (finished0==0)
  {
    // adapt the interface
    if (level<0)
      adaptp( *interface , metrics , params , interface_out , level+1 , finished );
    else
    {
      // don't adapt, just copy (for now)
      interface->points().copy(interface_out.points());
      interface_out.TopologyBase::copy(*interface);
    }

    if (rank == 0 )
    {
      printf("topology_out volume = %g\n",topology_out.volume());
      printf("interface volume = %g\n",interface->volume());
      printf("interface out volume = %g\n",interface_out.volume());

      #if 0
      for (std::map<index_t,index_t>::iterator it=interface2output.begin();it!=interface2output.end();++it)
        printf("interface point %lu maps to %lu\n",it->first,it->second);
      #endif
    }

    std::vector<index_t> interface_points( interface_out.points().nb() );
    for (index_t k=0;k<interface_out.points().nb();k++)
    {
      if (interface2output.find(k)!=interface2output.end())
        interface_points[k] = interface2output.at(k);
      else
      {
        interface_points[k] = topology_out.points().nb();
        index_t p = topology_out.points().nb();
        topology_out.points().create( interface_out.points()[k] );
        topology_out.points().set_entity( p , interface_out.points().entity(k) );
        topology_out.points().set_param( p , interface_out.points().u(k) );
        topology_out.points().set_fixed( p , interface_out.points().fixed(k) );
      }
    }

    // accumulate the adapted interface into the global output topology
    index_t offset = topology_out.points().nb();
    for (index_t k=0;k<interface_out.nb();k++)
    {
      std::vector<index_t> s = interface_out.get(k);
      for (index_t i=0;i<s.size();i++)
        s[i] = interface_points[s[i]];
      topology_out.add( s.data() , s.size() );
    }


    //topology_out.Table<index_t>::print();

    for (index_t k=0;k<interface_out.nb_children();k++)
      topology_out.add_child( interface_out.child_smptr(k) );
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

#endif
