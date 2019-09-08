// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXD_NORM_H
#define MATRIXD_NORM_H

#include "tools/SANSnumerics.h"
#include "tools/SANSException.h"
#include "MatrixD.h"

namespace SANS
{
namespace DLA
{

inline Real
normFrobenius(const MatrixDView<Real>& mat)
{
  Real sum = 0.0;

  for (int i = 0; i < mat.m(); i++)
    for (int j = 0; j < mat.n(); j++)
      sum += mat(i,j)*mat(i,j);

  return sqrt(sum);
}

} //namespace DLA
} //namespace SANS


#endif //MATRIXD_NORM_H
