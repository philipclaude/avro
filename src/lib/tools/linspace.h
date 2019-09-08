// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef LINSPACE_H
#define LINSPACE_H

#include "tools/SANSnumerics.h"     // Real
#include "tools/SANSException.h"

#include <vector>

namespace SANS
{

//returns a vector of consecutive integers in [a, b]
inline void
linspace( const int a, const int b, std::vector<int>& v )
{
  SANS_ASSERT_MSG( b >= a, "b = %d, a = %d", b, a);

  const int N = b - a + 1;

  v.resize(N);
  for (int i = 0; i < N; i++)
    v[i] = a + i;

}

inline std::vector<int>
linspace( const int a, const int b )
{
  std::vector<int> v;
  linspace( a, b, v );
  return v;
}

// Equivalent to MATLAB's linspace - generates a vector of equally spaced numbers between a and b.
inline std::vector<Real>
linspace( const Real a, const Real b, const int n )
{
  std::vector<Real> v(n);

  const Real L = b - a;

  for (int i = 0; i < n; i++)
    v[i] = a + ((Real)i / (Real)(n-1))*L;

  return v;
}

}
#endif // LINSPACE_H
