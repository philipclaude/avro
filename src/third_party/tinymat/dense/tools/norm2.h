// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef dense_NORM2_H
#define dense_NORM2_H

//based on LAPACK DLAPY2
//Compute sqrt(x*x + y*y), taking care not to cause unnecessary overflow
template< class T >
T norm2(const T& x, const T&y)
{
  T xabs = fabs(x);
  T yabs = fabs(y);
  T w = max(xabs,yabs);
  T z = min(xabs,yabs);
  if ( z == 0 )
    return w;
  else
    return w*sqrt( T(1) + (z/w)*(z/w) );
}

#endif //dense_NORM2_H
