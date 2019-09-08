// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define ILU0_INSTANTIATE
#include "ILU0_impl.h"
#include "LinearAlgebra/SparseLinAlg/SparseMatrix_CRS.h"

namespace SANS
{
namespace SLA
{


template class ILU0< SparseMatrix_CRS<Real> >;


} //namespace SLA
} //namespace SANS
