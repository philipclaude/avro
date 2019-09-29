// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATMUL_NATIVE_H
#define MATMUL_NATIVE_H

#include "tools/SANSnumerics.h"     // Real
#include "tools/minmax.h"

#include "tools/CacheLineSize.h"
#include "numpack/dense/dynamic/MatrixD_Type.h"

namespace numpack 
{
namespace DLA
{


template<class TL, class TR, class T>
class MatMul_Native
{
  static const int CacheItemsL = MAX(1, CACHE_LINE_SIZE / sizeof(TL) );
  static const int CacheItemsR = MAX(1, CACHE_LINE_SIZE / sizeof(TR) );

public:
//-----------------------------------------------------------------------------
  static void value(const MatrixDView<TL>& ML, const MatrixDView<TR>& MR, const Real sgn, MatrixDView<T>& res );

//-----------------------------------------------------------------------------
  static void plus(const MatrixDView<TL>& ML, const MatrixDView<TR>& MR, const Real sgn, MatrixDView<T>& res );

};


} //namespace DLA
} //namespace numpack 

#endif //MATMUL_NATIVE_H
