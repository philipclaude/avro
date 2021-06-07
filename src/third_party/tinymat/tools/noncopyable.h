// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef NONCOPYABLE_H
#define NONCOPYABLE_H

namespace tinymat 
{

// This is a convenient base class that can be used to prevent a class from being copied.
// See boost::noncopyable as well

class noncopyable
{
protected:
  noncopyable() {}
  ~noncopyable() {}

  noncopyable( const noncopyable& ) = delete;
  noncopyable operator=( const noncopyable& ) = delete;
};

}

#endif //NONCOPYABLE_H
