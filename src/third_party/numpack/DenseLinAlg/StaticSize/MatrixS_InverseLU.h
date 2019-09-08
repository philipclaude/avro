// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXS_INVERSELU_H
#define MATRIXS_INVERSELU_H

#include "tools/SANSnumerics.h"     // Real
#include "MatrixS_Type.h"
#include "MatrixS_Inverse.h"
#include "MatrixS_Decompose_LU.h"

//Specialization of matrix inverse to solve a linear system using LU decomposition without pivoting

namespace numpack 
{
namespace DLA
{

// Structure for solving using LU without pivoting
  template< int M, int N, class T, class MatrixType >
  struct MatrixSLUSolver
  {
    static_assert(M == N, "LU only works for square matrices");

    //Converts the matrix to an LU decomposed matrix without pivoting and unit diagonal on U
    typedef MatrixSLU<M, T> FactorType;

    static void Solve( const FactorType& Factorized, const Real sgn, MatrixType& res );
  };

} //namespace DLA
} //namespace numpack 

#endif //MATRIXS_INVERSELU_H
