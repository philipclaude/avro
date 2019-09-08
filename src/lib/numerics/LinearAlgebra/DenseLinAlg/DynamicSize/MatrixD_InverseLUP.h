// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXD_INVERSELUP_H
#define MATRIXD_INVERSELUP_H

#include "tools/SANSnumerics.h"     // Real
#include "MatrixD_Type.h"
#include "MatrixD_Inverse.h"
#include "MatrixD_Decompose_LUP.h"

//Specialization of matrix inverse to solve a linear system using LU decomposition with pivoting

namespace SANS
{
namespace DLA
{

// Structure for solving using LU with pivoting
template< class T, class MatrixType >
struct MatrixDLUPSolver
{
  //Convert the matrix to an LU decomposed matrix with pivoting
  typedef MatrixDLUP<T> FactorType;

  static void Solve( const FactorType& Matrix, const Real sgn, MatrixType& res );
};

} //namespace DLA
} //namespace SANS

#endif //MATRIXD_INVERSELUP_H
