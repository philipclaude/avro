// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXD_INVERSELU_H
#define MATRIXD_INVERSELU_H

#include "tools/SANSnumerics.h"     // Real
#include "MatrixD_Type.h"
#include "MatrixD_Inverse.h"
#include "MatrixD_Decompose_LU.h"

//Specialization of matrix inverse to solve a linear system using LU decomposition without pivoting

namespace tinymat 
{
namespace DLA
{

// Structure for solving using LU without pivoting
template< class T, class MatrixType >
struct MatrixDLUSolver
{
  //Converts the matrix to an LU decomposed matrix without pivoting and unit diagonal on U
  typedef MatrixDLU<T> FactorType;

  static void Solve( const FactorType& Factorized, const Real sgn, MatrixType& res );
};

} //namespace DLA
} //namespace tinymat 

#endif //MATRIXD_INVERSELU_H
