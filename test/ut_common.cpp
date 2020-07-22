#include "common/process.h"
#include "common/tools.h"

#include "numerics/predicates.h"

typedef avro::real_t REAL;
#include <tetgen1.5.0/predicates.h>

#include <triangle/predicates.h>

using namespace avro;

void
ut_pre(int argc, char** argv)
{
  UNUSED(argc);
  UNUSED(argv);

  print_info();

  // initialize the distributed/shared memory task/thread managers
  #ifdef AVRO_MPI
  ProcessMPI::initialize();
  #endif
  ProcessCPU::initialize();
  ProcessGPU::initialize();

  // initialize the predicates
  exactinit(1,0,0,10,10,10);
  exactinit();
  GEO::PCK::initialize();
}

void
ut_post()
{
  // nothing to do
}
