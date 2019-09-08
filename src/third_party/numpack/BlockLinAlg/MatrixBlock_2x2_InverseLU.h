// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXBLOCK_2X2_INVERSELU_H
#define MATRIXBLOCK_2X2_INVERSELU_H

#include "tools/SANSnumerics.h"     // Real

#include "MatrixBlock_2x2.h"
#include "VectorBlock_2.h"

#include "MatrixBlock_2x2_Inverse.h"
#include "MatrixBlock_2x2_Decompose_LU.h"

//Specialization of matrix inverse to solve a linear system using LU decomposition without pivoting

namespace numpack 
{
namespace BLA
{

// Structure for solving using LU without pivoting
  template< class M00, class M01, class M10, class M11, class MatrixType >
  struct MatrixBlock_2x2_LUSolver
  {
    //Converts the matrix to an LU decomposed matrix without pivoting and unit diagonal on U
    typedef MatrixBlock_2x2_LU<M00, M01, M10, M11> FactorType;

    static void Solve( const FactorType& Factorized, const Real sgn, MatrixType& res );
  };

} //namespace BLA
} //namespace numpack 

#endif //MATRIXBLOCK_2X2_INVERSELU_H
