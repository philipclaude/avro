// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXS_EXP_H
#define MATRIXS_EXP_H

#include <cmath>  // exp

#include "tools/SANSnumerics.h"

#include "MatrixS_Type.h"
#include "MatrixS_Transpose.h"
#include "MatrixS_Diag.h"
#include "Eigen.h"

//=============================================================================
template< int M, class T >
inline SANS::DLA::MatrixSymS<M,T>
exp(const SANS::DLA::MatrixSymS<M,T>& A)
{
  //Compute the eigensystem with normalized eigenvectors
  SANS::DLA::EigenSystemPair<M,T> LE(A);

  //Compute the exponential of the eigenvalues
  for (int i = 0; i < M; i++ )
    LE.L[i] = exp(LE.L[i]);

  //Return the Symmetric matrix
  return LE;
}

//=============================================================================
template< int M, class T >
inline SANS::DLA::EigenSystemPair<M,T>
exp(const SANS::DLA::EigenSystemPair<M,T>& LE)
{
  //Copy the Eigen pair
  SANS::DLA::EigenSystemPair<M,T> expLE(LE);

  //Compute the exponential of the eigenvalues
  for (int i = 0; i < M; i++ )
    expLE.L[i] = exp(LE.L[i]);

  //Return the Eigen pair
  return expLE;
}

//=============================================================================
template< int M, class T >
inline SANS::DLA::VectorS<M,T>
exp(const SANS::DLA::VectorS<M,T>& v)
{
  SANS::DLA::VectorS<M,T> vr;

  //Compute the exponential of each component
  for (int i = 0; i < M; i++ )
    vr[i] = exp(v[i]);

  return vr;
}

#endif //MATRIXS_EXP_H
