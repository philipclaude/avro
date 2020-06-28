#include "unit_tester.hpp"

#include "adaptation/parallel.h"

#include "common/process.h"

#include "geometry/egads/context.h"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/egads.h"
#include "library/plots.h"

#include "mesh/partition.h"
#include "mesh/points.h"
#include "mesh/topology.h"

using namespace avro;

UT_TEST_SUITE( mesh_mpi_suite )

#ifdef AVRO_MPI

UT_TEST_CASE( test1 )
{
  // this test case mostly just plays ping-pong with the sender/receiver of the mesh partitions
  coord_t number = 2;
  coord_t dim = number;
  coord_t udim = dim -1;

  mpi::communicator& comm = ProcessMPI::get_comm();
  int rank = mpi::rank();

  std::vector<index_t> dims(number,10);
  CKF_Triangulation topology(dims);
  topology.neighbours().fromscratch() = true;
  topology.neighbours().compute();

  EGADS::Context context;
  std::vector<real_t> lens(number,1.);
  EGADS::Cube geometry(&context,lens);
  topology.points().attach(geometry);

  std::vector<Entity*> entities;
  geometry.get_entities(entities);

  index_t nb_partition = TEST_NUM_PROCS -1;

  // interface topology
  AdaptationInterface<Simplex> interface(dim,udim,number);
  topology.points().copy( interface.points() );

  // send the partitions to each processor
  Topology_Partition<Simplex> topology_p(dim,udim,number);
  topology_p.set_entities(entities);
  if (rank == 0 )
  {
    // partition the topology
    Partition<Simplex> partition(topology);
    partition.compute(nb_partition);

    std::vector<std::set<index_t>> halo_points( nb_partition );
    partition.compute_interface_points(  halo_points );

    std::vector< std::shared_ptr<Topology_Partition<Simplex> > > parts(nb_partition);
    partition.get(parts);

    #if 0
    for (index_t k=0;k<parts.size();k++)
    {
      parts[k]->compute_crust();

      const std::vector<index_t>& crust = parts[k]->crust();

      // add the crust to the interface
      // everything is in local partition coordinates
      index_t np = interface.points().nb();
      printf("partition %lu, offset = %lu\n",k,np);
      for (index_t j=0;j<crust.size();j++)
      {
        //if (k>1) break;
        std::vector<index_t> s = parts[k]->get( crust[j] );
        print_inline(s,"interface before map: ");

        for (index_t i=0;i<s.size();i++)
          s[i] = parts[k]->local2global(s[i]);

        print_inline(s,"interface elem: ");
        interface.Topology<index_t>::add(s.data(),s.size());
      }

      parts[k]->move_to_front( parts[k]->halo() );
      parts[k]->remove_elements( parts[k]->crust() );
      parts[k]->remove_unused();
    }
    #else
    for (index_t k=0;k<parts.size();k++)
    {
      parts[k]->set_entities(entities);

      // extract the halo for this partition and map to local indices
      std::vector<index_t> halo( halo_points[k].begin() , halo_points[k].end() );
      for (index_t j=0;j<halo.size();j++)
        halo[j] = parts[k]->global2local( halo[j] );

      print_inline(halo);

      // compute the set of elements touching a halo point
      parts[k]->neighbours().fromscratch();
      parts[k]->neighbours().compute();
      parts[k]->inverse().build();
      std::vector<index_t> ball;
      std::vector<index_t> crust;
      for (index_t j=0;j<halo.size();j++)
      {
        ball.clear();
        parts[k]->inverse().ball( halo[j] , ball );
        for (index_t i=0;i<ball.size();i++)
          crust.push_back( ball[i] );
      }
      uniquify(crust);

      // add the crust to the interface
      // crust elements are in local indices, convert to global
      for (index_t j=0;j<crust.size();j++)
      {
        std::vector<index_t> s = parts[k]->get( crust[j] );
        for (index_t i=0;i<s.size();i++)
          s[i] = parts[k]->local2global(s[i]);
        interface.Topology<index_t>::add(s.data(),s.size());
      }

      printf("partition %lu has %lu elements with %lu elements in the crust\n",k+1,parts[k]->nb(),crust.size());

      // remove crust elements from the partition since they do not get adapted
      parts[k]->move_to_front( halo );
      parts[k]->remove_elements( crust );
      parts[k]->remove_unused();
    }
    #endif

    interface.remove_unused();

    printf("sending partitions..\n");
    for (index_t k=0;k<parts.size();k++)
    {
      UT_ASSERT( parts[k]->all_points_accounted() );
      parts[k]->inverse().build();
      parts[k]->send( comm , k+1 );
    }
  }
  else
  {
    printf("receiving partition..\n");
    avro_assert( topology_p.points().nb() == 0 );
    topology_p.receive( comm , 0 );
    UT_ASSERT( topology_p.all_points_accounted() );
    topology_p.inverse().build();
  }
  mpi::barrier();

  Points points(dim,udim);
  Topology<Simplex> topology_out(points,number);
  if (rank == 0)
  {
    for (index_t i=0;i<nb_partition;i++)
    {
      Topology_Partition<Simplex> tk(dim,udim,number);
      tk.set_entities(entities);
      tk.receive(comm,i+1);

      //tk.compute_mantle();

      index_t np = points.nb();
      for (index_t k=0;k<tk.points().nb();k++)
      {
        points.create( tk.points()[k] );
        points.set_entity( np+k , tk.points().entity(k) );
        for (coord_t d=0;d<udim;d++)
          points.u( np+k  )[d] = tk.points().u(k)[d];
      }

      for (index_t k=0;k<tk.nb();k++)
      {
        for (index_t j=0;j<tk.nv(k);j++)
          tk(k,j) += np;
        topology_out.add( tk(k) , tk.nv(k) );
      }
    }
  }
  else
  {
    topology_p.send( comm , 0 );
  }
  mpi::barrier();

  if (rank==0)
  {
    graphics::Visualizer vis;
    vis.add_topology(topology_out);
    vis.add_topology(interface);
    library::Plot<Simplex> plot(interface.points());
    vis.add_topology(plot);
    vis.run();
  }
  else
  {
    UT_ASSERT( points.nb() == 0 );
    UT_ASSERT( topology_out.nb() ==0 );
  }

}
UT_TEST_CASE_END( test1 )

#endif


UT_TEST_SUITE_END( mesh_mpi_suite )
