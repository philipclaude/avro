//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
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
  GEO::PCK::initialize();*/
}

void
ut_post()
{
  // nothing to do
}
