// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef CHKDENSEMATRIX_H
#define CHKDENSEMATRIX_H

#include <iostream>

//----------------------------------------------------------------------------//
template< class MatrixType >
bool chkMatrixD22( const MatrixType& z, int m11, int m12, int m21, int m22 )
{
  bool isEqual = true;
  if (z(0,0) != m11) isEqual = false;
  if (z(0,1) != m12) isEqual = false;
  if (z(1,0) != m21) isEqual = false;
  if (z(1,1) != m22) isEqual = false;
  if (!isEqual)
  {
    std::cout << "actual ("
              << "(" << z(0,0) << " " << z(0,1) << ") "
              << "(" << z(1,0) << " " << z(1,1) << ") "
              << ")" << std::endl;
    std::cout << "  expected ("
              << "(" << m11 << " " << m12 << ") "
              << "(" << m21 << " " << m22 << ") "
              << ")" << std::endl;
  }
  return isEqual;
}

//----------------------------------------------------------------------------//
template< class MatrixType >
bool chkMatrixD21( const MatrixType& z, int m11, int m21 )
{
  bool isEqual = true;
  if (z(0,0) != m11) isEqual = false;
  if (z(1,0) != m21) isEqual = false;
  if (!isEqual)
  {
    std::cout << "actual ("
              << "(" << z(0,0) << ") " << std::endl
              << "(" << z(1,0) << ") "
              << ")" << std::endl;
    std::cout << "  expected ("
              << "(" << m11 << ") " << std::endl
              << "(" << m21 << ") "
              << ")" << std::endl;
  }
  return isEqual;
}

//----------------------------------------------------------------------------//
template< class MatrixType >
bool chkMatrixD12( const MatrixType& z, int m11, int m12 )
{
  bool isEqual = true;
  if (z(0,0) != m11) isEqual = false;
  if (z(0,1) != m12) isEqual = false;
  if (!isEqual)
  {
    std::cout << "actual ("
              << "(" << z(0,0) << ") " << std::endl
              << "(" << z(0,1) << ") "
              << ")" << std::endl;
    std::cout << "  expected ("
              << "(" << m11 << ") " << std::endl
              << "(" << m12 << ") "
              << ")" << std::endl;
  }
  return isEqual;
}

//----------------------------------------------------------------------------//
template< class T, class Int>
bool chkMatrixSymD22( const numpack::DLA::MatrixSymD<T>& z, Int m00, Int m01, Int m10, Int m11 )
{
  if (z.m()!=2 || z.n()!=2) return false;

  bool isEqual = true;
  if (z(0,0) != m00) isEqual = false;
  if (z(0,1) != m01) isEqual = false;
  if (z(1,0) != m10) isEqual = false;
  if (z(1,1) != m11) isEqual = false;
  if (!isEqual)
  {
    std::cout << "actual   ("
              << "(" << z(0,0) << " " << z(0,1) << ") "
              << "(" << z(1,0) << " " << z(1,1) << ") "
              << ")" << std::endl;
    std::cout << "expected ("
              << "(" << m00 << " " << m01 << ") "
              << "(" << m10 << " " << m11 << ") "
              << ")" << std::endl;
  }
  return isEqual;
}

//----------------------------------------------------------------------------//
template< class T, class Int>
bool chkVectorD2( const numpack::DLA::VectorD<T>& z, Int a, Int b )
{
  if (z.m()!=2){
    printf("bad size!!\n");
    return false;
  }

  bool isEqual = true;
  if ((z[0] != a) || (z[1] != b))
  {
    isEqual = false;
    std::cout << "actual (" << z << ")  expected ("
              << "(" << a << " " << b << "))" << std::endl;
  }
  return isEqual;
}

#endif //CHKDENSEMATRIX_H
