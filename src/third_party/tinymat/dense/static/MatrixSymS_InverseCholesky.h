// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXSYMS_INVERSECHOLESKY_H
#define MATRIXSYMS_INVERSECHOLESKY_H

#include "tools/SANSnumerics.h"     // Real
#include "MatrixS_Type.h"
#include "MatrixS_Inverse.h"
#include "MatrixSymS_Decompose_Cholesky.h"

//Specialization of matrix inverse to solve a linear system using Cholesky decomposition

namespace tinymat 
{
namespace DLA
{

  // Structure for solving using Cholesky on a symmetric matrix
  template< int M, int N, class T, class MatrixType >
  struct MatrixSymSCholeskySolver
  {
    static_assert(M == N, "Cholesky only works for square matrices");

    //Converts the symmetric positive definite matrix to an Cholesky decomposition L*L^T
    typedef MatrixSymSCholesky<M, T> FactorType;

    static void Solve( const FactorType& Factorized, const Real sgn, MatrixType& res );
  };

} //namespace DLA
} //namespace tinymat 

#endif //MATRIXSYMS_INVERSECHOLESKY_H
