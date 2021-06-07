// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef TOOLS_REF_H
#define TOOLS_REF_H

namespace tinymat 
{

// A simple class to store a reference to an object

template<class T>
struct Ref
{
public:
  // cppcheck-suppress noExplicitConstructor
  Ref(const T& r) : ref(r) {}

  operator const T&() const { return ref; }

  const T& ref;
};

}

#endif //TOOLS_REF_H
