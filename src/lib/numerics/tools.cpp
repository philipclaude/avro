// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "numerics/expansion.h"
#include "numerics/predicates.h"
#include "numerics/tools.h"

#include <cmath>

namespace avro
{

namespace numerics
{

real_t
sum( const std::vector<real_t>& x ) {
  // implements the Kahan sum of the values in a vector
  real_t sum = 0.;
  real_t c = 0.;
  for (index_t k = 0; k < x.size(); k++) {
    real_t y = x[k] -c;
    real_t t = sum +y;
    c = (t -sum) -y;
    sum = t;
  }
  return sum;
}

real_t
exactsum( const std::vector<real_t>& x ) {
  using namespace GEO;
  using namespace GEO::PCK;
  real_t sum = 0.;
  for (index_t k = 0;k < x.size(); k++) {
    expansion& s = expansion_sum(sum,x[k]);
    sum = s.data()[s.length()-1];
    //printf("accumulated %.20e\n",sum);
  }
  return sum;
}

real_t
naivesum( const std::vector<real_t>& x ) {
  real_t sum = 0.;
  for (index_t k = 0; k < x.size(); k++)
    sum += x[k];
  return sum;
}

real_t
norm( const std::vector<real_t>& x ) {
  real_t sum = 0.;
  for (index_t k = 0; k < x.size(); k++)
    sum += x[k]*x[k];
  return std::sqrt(sum);
}

} // numerics

} // avro
