// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXSYMS_INVERSELDLT_H
#define MATRIXSYMS_INVERSELDLT_H

#include "tools/SANSnumerics.h"     // Real
#include "MatrixS_Type.h"
#include "MatrixS_Inverse.h"
#include "MatrixSymS_Decompose_LDLT.h"

//Specialization of matrix inverse to solve a linear system using LDL^TT decomposition

namespace SANS
{
namespace DLA
{

  // Structure for solving using LDLT on a symmetric matrix
  template< int M, int N, class T, class MatrixType >
  struct MatrixSymSLDLTSolver
  {
    static_assert(M == N, "LDL^T only works for square matrices");

    //Converts the symmetric matrix to an LDL^T decomposition
    typedef MatrixSymSLDLT<M, T> FactorType;

    static void Solve( const FactorType& Factorized, const Real sgn, MatrixType& res );
  };

} //namespace DLA
} //namespace SANS

#endif //MATRIXSYMS_INVERSELDLT_H
