// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define EIGEN_INSTANTIATE
#define LAPACK_GEEV LAPACK_dgeev

//Instantiate Eigen value routines for double precision

#ifdef DLA_LAPACK

#include "Eigen_impl.h"

namespace numpack 
{
namespace DLA
{
template struct LAPACK_Eigen<double>;
}
}

#endif