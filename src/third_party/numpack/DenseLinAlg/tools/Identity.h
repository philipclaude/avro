// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef DENSEMATRIX_IDENTITY_H
#define DENSEMATRIX_IDENTITY_H

namespace numpack 
{
namespace DLA
{

//This a datatype used to indicate that a matrix should be set to Identity
struct Identity
{
  operator int() const { return 1; }
};

} //namespace DLA
} //namespace numpack 


#endif //DENSEMATRIX_IDENTITY_H
