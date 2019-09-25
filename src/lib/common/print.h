// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef URSA_COMMON_PRINT_H_
#define URSA_COMMON_PRINT_H_

#include "common/types.h"

#include <vector>
#include <string>
#include <stdarg.h>

namespace ursa
{

void ursaPrintf( const char* fmt , ... );

template<typename type>
static void
printValue( const type& x );

template<>
inline void
printValue( const index_t& x )
{
    printf("%d ",int(x));
}

template<>
inline void
printValue( const int& x )
{
    printf("%d ",int(x));
}


template<>
inline void
printValue( const real& x )
{
    printf("%g ",x);
}

template<typename type>
static void
printInline( const std::vector<type>& s , const std::string& name=std::string() , const int id=-1 , const index_t nt=0 )
{
  for (index_t k=0;k<nt;k++) printf("  ");
  if (name.size()>0)
      printf("%s",name.c_str());
  if (id>=0)
      printf("[%d]: ",id);
  printf("( ");
  for (index_t j=0;j<s.size();j++)
      printValue(s[j]);
  printf(")\n");
}

inline void
tabbedPrint( const index_t nt , const char *fmt , ... )
{
  for (index_t k=0;k<nt;k++) printf("  ");
  va_list args;
  va_start(args,fmt);
  char buffer[256];
  vsprintf(buffer,fmt,args);
  va_end(args);
}

template void printInline( const std::vector<index_t>& s , const std::string& name , const int id , const index_t nt );
template void printInline( const std::vector<real>& s , const std::string& name , const int id , const index_t nt );

} // ursa

#endif
