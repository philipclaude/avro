// ursa: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef URSA_COMMON_INFO_H_
#define URSA_COMMON_INFO_H_

#include <cstdio>

inline void
printInfo()
{

  printf("\nursa compiled with ");

  if (__cplusplus == 201103L) printf("c++11");
  else if (__cplusplus == 199711L ) printf("c++98");
  else printf("pre-standard c++");

  printf("deleting memory upon termination using");
  #ifdef URSA_SMART_PTR
    printf(" internal smart pointers.\n");
  #else
    printf(" std::shared_ptr (c++11).\n");
  #endif
  printf("\n");
}

#endif
