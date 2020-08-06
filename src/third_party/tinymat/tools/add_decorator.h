// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef ADD_DECORATOR_H
#define ADD_DECORATOR_H

namespace tinymat 
{

// Used to apply a decorator class in template meta functions
template<template<class> class decorator, class arg>
struct add_decorator
{
  template<class arg1>
  struct apply
  {
    typedef decorator<arg1> type;
  };
};

// Used to apply an ND decorator class in template meta functions
template<template<class,class> class decorator, class arg>
struct add_ND_decorator
{
  template<class arg1>
  struct apply
  {
    typedef decorator<typename arg1::PhysDim, arg1> type;
  };
};

} //namespace tinymat 

#endif //ADD_DECORATOR_H
