// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define LU_SGS_INSTANTIATE
#include "LU_SGS_impl.h"
#include "numpack/sparse/SparseMatrix_CRS.h"
#include "numpack/dense/static/VectorS.h"

namespace numpack 
{
namespace SLA
{


template class LU_SGS< SparseMatrix_CRS< DLA::MatrixS<1,1,Real> > >;
template class LU_SGS< SparseMatrix_CRS< DLA::MatrixS<2,2,Real> > >;
template class LU_SGS< SparseMatrix_CRS< DLA::MatrixS<4,4,Real> > >;


} //namespace SLA
} //namespace numpack 
