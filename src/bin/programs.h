#ifndef AVRO_BIN_PROGRAMS_H_
#define AVRO_BIN_PROGRAMS_H_

#include "common/tools.h"
#include "avro_types.h"

#include <memory>
#include <stdlib.h>
#include <string>
#include <vector>

namespace avro
{

namespace numerics
{
class DiscreteField;
}

class Model;
class Points;
class Mesh;
template<typename type> class Topology;

namespace programs
{

void help();
std::string lookfor( const char** args , int nb_args , const std::string& option );

template<typename type>
bool
parse( const std::string& str_arg , type& arg )
{
  if (str_arg.empty()) return false; // do not modify the argument
  arg = unstringify<type>(str_arg);
  return true;
}

int adapt( int argc , const char** argv );
int plot( int argc , const char** argv , bool webplot=false );
int voronoi( int argc , const char** argv );
int conformity( int argc , const char** argv );
int convert( int argc , const char** argv );

} // programs

} // avro

#endif
