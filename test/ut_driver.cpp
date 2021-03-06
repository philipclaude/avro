//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifdef STANDALONE
#error "this should be used for entire unit tests"
#endif

#include "unit_tester.hpp"

#include "common/tools.h"
#include "common/process.h"

#include "numerics/predicates.h"

typedef avro::real_t REAL;
#include <tetgen1.5.0/predicates.h>

#include <triangle/predicates.h>

TestDriver* __driver__ = new TestDriver;

using namespace avro;

int
main()
{
  initialize_avro();

  // initialize the distributed/shared memory task/thread managers
  /*
  #if AVRO_MPI
  ProcessMPI::initialize();
  #endif
  ProcessCPU::initialize();
  ProcessGPU::initialize();

  // initialize the predicates
  exactinit(1,0,0,10,10,10);
  exactinit();
  GEO::PCK::initialize();
  */

  #ifdef STDOUT_REDIRECT
  FILE *fid = fopen(STDOUT_REDIRECT,"w");
  fclose(fid);
  #endif
  TestResult result;
  __driver__->run(result);
  result.summary();
  delete __driver__;
  if (!result.successful()) return 1;
  return 0;
}
