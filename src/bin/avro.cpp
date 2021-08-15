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

#include <string>

using namespace avro;

int
main( int argc , const char** argv )
{
  // check the global calling sequence
  if (argc<2)
  {
    printf("usage:\n\tavro -[program] [inputs]\n");
    programs::help();
    return 1;
  }

  // check the hyphen delimiter exists
  std::string hyphen(argv[1],argv[1]+1);
  if (hyphen!="-")
  {
    printf("forgot a hyphen in front of program name!\n");
    printf("usage:\n\tavro -[program] [inputs]\n");
    programs::help();
    return 1;
  }

  // get the name of the program
  std::string program(argv[1]+1);

  // only pass in relevant arguments to the program
  int nb_inputs = argc -2;
  const char** inputs = argv +2;

  printf("running program %s\n",program.c_str());
  int result = 0;
  if (program=="adapt")
    result = programs::adapt(nb_inputs,inputs);
  else if (program=="plot")
    result = programs::plot(nb_inputs,inputs);
  else if (program=="webplot")
    result = programs::plot(nb_inputs,inputs,true);
  else if (program=="conformity")
    result = programs::conformity(nb_inputs,inputs);
  else if (program=="convert")
    result = programs::convert(nb_inputs,inputs);
  else if (program=="check")
    result = programs::check(nb_inputs,inputs);
  else if (program == "animate")
    result = programs::animate(nb_inputs,inputs);
  else {
    printf("unknown program :(\n");
    result = 1;
  }

  return result;
}
