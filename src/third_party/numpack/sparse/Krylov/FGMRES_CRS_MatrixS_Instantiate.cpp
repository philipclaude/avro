// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define FGMRES_INSTANTIATE
#include "FGMRES_impl.h"

#include "numpack/dense/static/VectorS.h"

namespace numpack 
{
namespace SLA
{

template class FGMRES< SparseMatrix_CRS< DLA::MatrixS<1,1,Real> > >;
template class FGMRES< SparseMatrix_CRS< DLA::MatrixS<2,2,Real> > >;
template class FGMRES< SparseMatrix_CRS< DLA::MatrixS<4,4,Real> > >;
template class FGMRES< SparseMatrix_CRS< DLA::MatrixS<6,6,Real> > >;

}
}