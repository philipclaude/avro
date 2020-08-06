// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXS_POW_H
#define MATRIXS_POW_H

#include <cmath> // pow

#include "tools/SANSnumerics.h"

#include "MatrixS_Type.h"
#include "MatrixS_Transpose.h"
#include "MatrixS_Diag.h"
#include "Eigen.h"

//=============================================================================
template< int M, class T, class Texp >
inline tinymat::DLA::MatrixSymS<M,T>
pow(const tinymat::DLA::MatrixSymS<M,T>& A, Texp r)
{
  //Compute the eigensystem with normalized eigenvectors
  tinymat::DLA::EigenSystemPair<M,T> LE(A);

  //Compute the power of the eigenvalues
  for (int i = 0; i < M; i++ )
    LE.L[i] = pow(LE.L[i], r);

  //Return the Symmetric matrix
  return LE;
}

//=============================================================================
template< int M, class T, class Texp >
inline tinymat::DLA::EigenSystemPair<M,T>
pow(const tinymat::DLA::EigenSystemPair<M,T>& LE, Texp r)
{
  //Copy the Eigen pair
  tinymat::DLA::EigenSystemPair<M,T> powLE(LE);

  //Compute the power of the eigenvalues
  for (int i = 0; i < M; i++ )
    powLE.L[i] = pow(LE.L[i], r);

  //Return the Eigen pair
  return powLE;
}

//=============================================================================
template< int M, class T >
inline tinymat::DLA::MatrixSymS<M,T>
sqrt(const tinymat::DLA::MatrixSymS<M,T>& A)
{
  //Compute the eigensystem with normalized eigenvectors
  tinymat::DLA::EigenSystemPair<M,T> LE(A);

  //Compute the sqrt of the eigenvalues
  for (int i = 0; i < M; i++ )
    LE.L[i] = sqrt(LE.L[i]);

  //Return the Symmetric matrix
  return LE;
}

//=============================================================================
template< int M, class T >
inline tinymat::DLA::EigenSystemPair<M,T>
sqrt(const tinymat::DLA::EigenSystemPair<M,T>& LE)
{
  //Copy the Eigen pair
  tinymat::DLA::EigenSystemPair<M,T> sqrtLE(LE);

  //Compute the sqrt of the eigenvalues
  for (int i = 0; i < M; i++ )
    sqrtLE.L[i] = sqrt(LE.L[i]);

  //Return the Eigen pair
  return sqrtLE;
}

//=============================================================================
template< int M, class T, class Texp >
inline tinymat::DLA::VectorS<M,T>
pow(const tinymat::DLA::VectorS<M,T>& v, Texp r)
{
  tinymat::DLA::VectorS<M,T> vr;

  //Compute the power of each component
  for (int i = 0; i < M; i++ )
    vr[i] = pow(v[i], r);

  return vr;
}

#endif //MATRIXS_POW_H
