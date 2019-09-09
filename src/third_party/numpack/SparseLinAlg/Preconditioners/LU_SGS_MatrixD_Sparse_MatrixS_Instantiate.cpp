// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define LU_SGS_INSTANTIATE
#include "LU_SGS_impl.h"
#include "numpack/SparseLinAlg/SparseMatrix_CRS.h"
#include "numpack/DenseLinAlg/StaticSize/VectorS.h"

namespace numpack 
{
namespace SLA
{


template class LU_SGS< DLA::MatrixD< SparseMatrix_CRS< Real > > >;
template class LU_SGS< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixS<1,1,Real> > > >;
template class LU_SGS< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixS<2,2,Real> > > >;
template class LU_SGS< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixS<3,3,Real> > > >;
template class LU_SGS< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixS<4,4,Real> > > >;
template class LU_SGS< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixS<5,5,Real> > > >;
template class LU_SGS< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixS<6,6,Real> > > >;
template class LU_SGS< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixS<7,7,Real> > > >;
template class LU_SGS< DLA::MatrixD< SparseMatrix_CRS< DLA::MatrixS<8,8,Real> > > >;


} //namespace SLA
} //namespace numpack 