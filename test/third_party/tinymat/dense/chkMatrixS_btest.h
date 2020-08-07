// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef CHKMATRIXS_BTEST_H
#define CHKMATRIXS_BTEST_H

#include "tinymat/dense/static/MatrixS.h"
#include "tinymat/dense/static/MatrixSymS.h"
#include "tinymat/dense/static/VectorS.h"

//----------------------------------------------------------------------------//
template< class T, class Int>
bool chkMatrixS22( const tinymat::DLA::MatrixS<2,2,T>& z, Int m00, Int m01, Int m10, Int m11 )
{
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
bool chkMatrixSymS22( const tinymat::DLA::MatrixSymS<2,T>& z, Int m00, Int m01, Int m10, Int m11 )
{
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
bool chkVectorS2( const tinymat::DLA::VectorS<2,T>& z, Int a, Int b )
{
  bool isEqual = true;
  if ((z[0] != a) || (z[1] != b))
  {
    isEqual = false;
    std::cout << "actual (" << z << ")  expected ("
              << "(" << a << " " << b << "))" << std::endl;
  }
  return isEqual;
}


#endif //CHKMATRIXS_BTEST_H
