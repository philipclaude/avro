#ifdef STANDALONE
#error "this should be used for entire unit tests"
#endif

#include "unit_tester.hpp"

#include "common/tools.h"
#include "common/process.h"

typedef avro::real_t REAL;
#include <tetgen1.5.0/predicates.h>

#include <triangle/predicates.h>

TestDriver* __driver__ = new TestDriver;

using namespace avro;

int
main()
{
  printInfo();

  // initialize the distributed/shared memory task/thread managers
  #ifdef avro_WITH_MPI
  ProcessMPI::initialize();
  #endif
  ProcessCPU::initialize();
  ProcessGPU::initialize();

  // initialize the predicates
  exactinit(1,0,0,10,10,10);
  exactinit();

  // create a graphics plotter instance
  //Server server;
  //Plotter* plotter = new Plotter(&server);
  //__plotter__ = (void*) plotter;
  //plotter->on() = false;

  #ifdef STDOUT_REDIRECT
  FILE *fid = fopen(STDOUT_REDIRECT,"w");
  fclose(fid);
  #endif
  TestResult result;
  __driver__->run(result);
  result.summary();
  delete __driver__;
  return 0;
}
