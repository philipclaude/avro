// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXS_LOG_H
#define MATRIXS_LOG_H

#include <cmath> // log

#include "tools/SANSnumerics.h"

#include "MatrixD_Type.h"
#include "MatrixD_Transpose.h"
#include "MatrixD_Diag.h"
#include "Eigen.h"

//=============================================================================
template< class T >
inline numpack::DLA::MatrixSymD<T>
log(const numpack::DLA::MatrixSymD<T>& A)
{
  //Compute the eigensystem with normalized eigenvectors
  numpack::DLA::EigenSystemPair<T> LE(A);

  //Compute the log of the eigenvalues
  for (int i = 0; i < A.m(); i++ )
    LE.L[i] = log(LE.L[i]);

  //Return the Symmetric matrix
  return LE;
}

#endif //MATRIXS_LOG_H
