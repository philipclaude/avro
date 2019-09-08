// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXBLOCK_2X2_TRANSPOSE_H
#define MATRIXBLOCK_2X2_TRANSPOSE_H

#include "tools/SANSnumerics.h"
#include "tools/SANSException.h"
#include "BlockLinAlg_Type.h"

#include "numpack/Transpose.h"

namespace numpack 
{
namespace BLA
{

//Represents the transpose of a 2x2 block matrix
template<class M00, class M01, class M10, class M11>
class MatrixBlock_2x2_Transpose : public BlockLinAlgType< MatrixBlock_2x2_Transpose<M00, M01, M10, M11> >
{

public:
  //Types of the transposed blocks
  typedef typename TransposeViewTraits<M00>::type MT00;
  typedef typename TransposeViewTraits<M10>::type MT01;
  typedef typename TransposeViewTraits<M01>::type MT10;
  typedef typename TransposeViewTraits<M11>::type MT11;

  explicit MatrixBlock_2x2_Transpose( const MatrixBlock_2x2<M00, M01, M10, M11>& M )
  : M_(M) {}

  // Lazy expression operations
  template<class T00, class T01, class T10, class T11>
  inline void value(const Real sgn, MatrixBlock_2x2<T00, T01, T10, T11>& res) const
  {
    res.m00 = sgn*Transpose(M_.m00);
    res.m01 = sgn*Transpose(M_.m10);
    res.m10 = sgn*Transpose(M_.m01);
    res.m11 = sgn*Transpose(M_.m11);
  }

  template<class T00, class T01, class T10, class T11>
  inline void plus(const Real sgn, MatrixBlock_2x2<T00, T01, T10, T11>& res) const
  {
    res.m00 += sgn*Transpose(M_.m00);
    res.m01 += sgn*Transpose(M_.m10);
    res.m10 += sgn*Transpose(M_.m01);
    res.m11 += sgn*Transpose(M_.m11);
  }

  int m() const { return M_.n(); } //Transpose, so m = M_.n
  int n() const { return M_.m(); } //Transpose, so n = M_.m

private:
  const MatrixBlock_2x2<M00, M01, M10, M11>& M_; //The transposed matrix
};

} //namespace BLA
} //namespace numpack 



#endif //MATRIXBLOCK_2X2_TRANSPOSE_H
