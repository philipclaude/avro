#include "unit_tester.hpp"

#include "adaptation/parallel.h"
#include "adaptation/parameters.h"

#include "common/mpi.hpp"
#include "common/process.h"

#include "graphics/application.h"

#include "library/ckf.h"

using namespace avro;

UT_TEST_SUITE( adaptation_parallel_test_suite )

UT_TEST_CASE( test1 )
{
  coord_t number = 2;

  mpi::communicator& comm = ProcessMPI::get_comm();
  int rank = mpi::rank();

  printf("creating topology\n");
  std::vector<index_t> dims(number,10);
  CKF_Triangulation topology(dims);

  topology.neighbours().fromscratch() = true;
  topology.neighbours().compute();

  Points points_out;
  Topology<Simplex> topology_out(points_out,number);
  AdaptationParameters params;

  params.nb_partition() = TEST_NUM_PROCS;

  std::vector<VertexMetric> metrics;
  int result = adaptp( topology , metrics , params , topology_out );



}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( adaptation_parallel_test_suite )
