//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "programs.h"

#include "common/directory.h"
#include "common/process.h"
#include "common/tools.h"

#include "graphics/application.h"

#include "library/factory.h"
#include "library/library.h"

#include "mesh/mesh.h"

#include "numerics/predicates.h"

#include "voronoi/cvt.h"

#include "avro.h"

typedef avro::real_t REAL;
#include <tetgen1.5.0/predicates.h>

#include <triangle/predicates.h>

#include <stdio.h>

namespace avro
{

namespace programs
{

int
sdot( int nb_input , const char** inputs )
{

  // so far only simplex adaptation is supported
  typedef Polytope type;

  if (nb_input<4 || nb_input==-1)
  {
    printf("\t\tvoronoi [input mesh] [sites] [geometry] [optional]\n");
    printf("\t\t--> optional can be:\n");
    printf("\t\t\tnb_iter=(int, # CVT iterations)\n");
    printf("\t\t\tnb_sites=(int, # sites for random only)\n");
    return 1;
  }

  ProcessCPU::initialize();
  ProcessGPU::initialize();

  // initialize the predicates
  exactinit(1,0,0,10,10,10);
  exactinit();
  GEO::PCK::initialize();

  // options
  bool found; UNUSED(found);
  const char **options = inputs +3;
  int  nb_options = nb_input -3;

  // get the input metric
  std::string sitesname( inputs[1] );

  return 0;
}


} // program

} // avro
