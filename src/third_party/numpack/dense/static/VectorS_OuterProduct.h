// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef VECTORS_OUTERPRODUCT_H
#define VECTORS_OUTERPRODUCT_H

#include "MatrixS_Type.h"
#include "numpack/dense/tools/PromoteSurreal.h"

namespace numpack 
{
namespace DLA
{
  //Computes the outer product of two VectorS
  template<int D, class Tx, class Ty>
  inline MatrixS<D,D,typename promote_Surreal<Tx,Ty>::type>
  OuterProduct(const VectorS<D,Tx>& x, const VectorS<D,Ty>& y)
  {
    typedef typename promote_Surreal<Tx,Ty>::type T;
    MatrixS<D,D,T> M;

    for (int i = 0; i < D; i++)
      for (int j = 0; j < D; j++)
        M(i,j) = x(i)*y(j);

    return M;
  }

  // computes the outer product of a vector with itself
  template<int D, class Tx>
  inline MatrixSymS<D,Tx>
  OuterProduct(const VectorS<D,Tx>& x)
  {
    MatrixSymS<D,Tx> M;

    for (int i = 0; i < D; i++)
      for (int j = i; j < D; j++) // only need to compute upper diagonal
        M(i,j) = x(i)*x(j);

    return M;
  }


} //namespace DLA
} //namespace numpack 

#endif //VECTORS_OUTERPRODUCT_H
