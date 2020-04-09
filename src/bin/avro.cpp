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
  else
  {
    printf("unknown program :(\n");
    result = 1;
  }

  return result;
}
