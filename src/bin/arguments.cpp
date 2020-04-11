#include "programs.h"

#include "common/tools.h"

namespace avro
{

namespace programs
{

// arguments help
void
help()
{
  // call each program to print the call usage
  adapt(-1,0);
  plot(-1,0);
  voronoi(-1,0);
}

// look for an argument in the list and return the value
std::string
lookfor( const char** args , int nb_args , const std::string& option )
{
  for (int k=0;k<nb_args;k++)
  {
    std::string s(args[k]);

    std::vector<std::string> S = split(s,"=");
    if (S.size()!=2)
    {
      printf("option %s not understood\n",args[k]);
      continue;
    }

    if (S[0]==option) return S[1];
  }
  return std::string();
}

} // programs

} // avro
