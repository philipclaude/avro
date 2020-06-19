#include "unit_tester.hpp"

#include "common/process.h"

#include "geometry/egads/context.h"

#include "library/ckf.h"
#include "library/egads.h"

#include "mesh/partition.h"
#include "mesh/points.h"
#include "mesh/topology.h"

using namespace avro;

UT_TEST_SUITE( mesh_mpi_suite )

#ifdef AVRO_MPI

UT_TEST_CASE( test1 )
{
  coord_t number = 3;
  coord_t dim = number;
  coord_t udim = dim -1;

  mpi::communicator& comm = ProcessMPI::get_comm();
  int rank = mpi::rank();

  std::vector<index_t> dims(number,4);
  CKF_Triangulation topology(dims);
  topology.neighbours().fromscratch() = true;
  topology.neighbours().compute();

  EGADS::Context context;
  std::vector<real_t> lens(number,1.);
  EGADS::Cube geometry(&context,lens);
  topology.points().attach(geometry);

  index_t nb_partition = TEST_NUM_PROCS -1;

  // send the partitions to each processor
  Topology_Partition<Simplex> topology_p(dim,udim,number);
  if (rank == 0 )
  {
    // partition the topology
    Partition<Simplex> partition(topology);
    partition.compute(nb_partition);

    std::vector< std::shared_ptr<Topology_Partition<Simplex> > > parts(nb_partition);
    partition.get(parts);

    printf("sending partitions..\n");
    for (index_t k=0;k<parts.size();k++)
    {
      parts[k]->send( comm , k+1 );
    }
  }
  else
  {
    printf("receiving partition..\n");
    topology_p.receive( comm , 0 );
  }
  mpi::barrier();

  if (rank == 0)
  {

    Points points(dim,udim);
    Topology<Simplex> topology_out(points,number);
    for (index_t i=0;i<nb_partition;i++)
    {
      Topology_Partition<Simplex> tk(number,number-1,number);
      tk.receive(comm,i+1);

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

    topology_out.Table<index_t>::print();
    topology_out.points().print();
  }
  else
  {
    topology_p.send( comm , 0 );
  }
  mpi::barrier();

}
UT_TEST_CASE_END( test1 )

#endif


UT_TEST_SUITE_END( mesh_mpi_suite )
