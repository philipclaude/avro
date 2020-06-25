#include "unit_tester.hpp"

#include "adaptation/parallel.h"
#include "adaptation/parameters.h"

#include "common/mpi.hpp"
#include "common/process.h"

#include "geometry/egads/context.h"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/egads.h"

using namespace avro;

UT_TEST_SUITE( adaptation_parallel_test_suite )

UT_TEST_CASE( test1 )
{
  coord_t number = 2;

  printf("creating topology\n");
  std::vector<index_t> dims(number,10);
  CKF_Triangulation topology(dims);

  EGADS::Context context;
  std::vector<real_t> lens(number,1.);
  EGADS::Cube geometry(&context,lens);
  topology.points().attach(geometry);

  topology.neighbours().fromscratch() = true;
  topology.neighbours().compute();

  Points points_out;
  Topology<Simplex> topology_out(points_out,number);
  AdaptationParameters params;
  params.standard();

  params.nb_partition() = TEST_NUM_PROCS -1;
  params.partition_size() = 50;

  std::vector<VertexMetric> metrics;
  int result = adaptp( topology , metrics , params , topology_out );



}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( adaptation_parallel_test_suite )
