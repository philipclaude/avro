// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef VECTORS_INNERPRODUCT_H
#define VECTORS_INNERPRODUCT_H

#include "MatrixS_Type.h"
#include "tinymat/dense/tools/PromoteSurreal.h"

namespace tinymat 
{
namespace DLA
{
  //Computes the inner product of two VectorS with a MatrixS: Transpose(x)*A*y
  template<int M, class Tx, class Ty, class Ta>
  inline typename promote_Surreal<Tx,Ty,Ta>::type
  InnerProduct(const VectorS<M,Tx>& x, const VectorS<M,Ty>& y, const MatrixS<M, M, Ta>& A)
  {
    typedef typename promote_Surreal<Tx,Ty,Ta>::type T;
    T sum = 0.0;
    for (int i = 0; i < M; i++)
      for (int j = 0; j < M; j++)
        sum += x(i)*A(i,j)*y(j);

    return sum;
  }

  //Computes the weighted inner product of a VectorS with itself: Transpose(x)*A*x
  template<int M, class Tx, class Ta>
  inline typename promote_Surreal<Tx,Ta>::type
  InnerProduct(const VectorS<M,Tx>& x, const MatrixS<M, M, Ta>& A)
  {
    typedef typename promote_Surreal<Tx,Ta>::type T;
    T sum = 0.0;
    for (int i = 0; i < M; i++)
      for (int j = 0; j < M; j++)
        sum += x(i)*A(i,j)*x(j);

    return sum;
  }

  //Computes the inner product of two VectorS with a MatrixSymS: Transpose(x)*A*y
  template<int M, class Tx, class Ty, class Ta>
  inline typename promote_Surreal<Tx,Ty,Ta>::type
  InnerProduct(const VectorS<M,Tx>& x, const VectorS<M,Ty>& y, const MatrixSymS<M, Ta>& A)
  {
    typedef typename promote_Surreal<Tx,Ty,Ta>::type T;
    T sum = 0.0;
    for (int i = 0; i < M; i++)
      for (int j = 0; j < M; j++)
        sum += x(i)*A(i,j)*y(j);

    return sum;
  }

  // Compute the weighted inner product of a VectorS with itself: Transpose(x)*A*y
  template<int M, class Tx, class Ta>
  inline typename promote_Surreal<Tx,Ta>::type
  InnerProduct(const VectorS<M,Tx>& x, const MatrixSymS<M, Ta>& A)
  {
    typedef typename promote_Surreal<Tx,Ta>::type T;
    T sum = 0.0;

    // for (int i = 0; i < M; i++)
    //   for (int j = 0; j < M; j++)
    //     sum += x(i)*A(i,j)*y(j);

    // faster version of the above, exploits symmetry
    for (int i = 0; i < M; i++)
    {
      sum += x(i)*A(i,i)*x(i); // diagonal
      for (int j = 0; j < i; j++)
        sum += 2*x(i)*A(i,j)*x(j); // lower diagonal
    }
    return sum;
  }


  //Computes the inner product of two VectorS: Transpose(x)*y, identical to DLA::Dot
  template<int M,class Tx, class Ty>
  inline typename promote_Surreal<Tx,Ty>::type
  InnerProduct(const VectorS<M,Tx>& x, const VectorS<M,Ty>& y)
  {
    typedef typename promote_Surreal<Tx,Ty>::type T;
    T sum = 0.0;
    for (int i = 0; i < M; i++)
        sum += x(i)*y(i);

    return sum;
  }

  //Computes the inner product of two VectorS with a third VectorS representing the diagonal: Transpose(x)*Diag(d)*y
  template<int M, class Tx, class Ty, class Td>
  inline typename promote_Surreal<Tx,Ty,Td>::type
  InnerProductDiag(const VectorS<M,Tx>& x, const VectorS<M,Ty>& y, const VectorS<M,Td>& d)
  {
    typedef typename promote_Surreal<Tx,Ty,Td>::type T;
    T sum = 0.0;
    for (int i = 0; i < M; i++)
        sum += x(i)*d(i)*y(i);

    return sum;
  }

  //Computes the inner product of a VectorS with itself, weighted by a diagonal matrix: Transpose(x)*Diag(d)*x
  template<int M, class Tx, class Td>
  inline typename promote_Surreal<Tx,Td>::type
  InnerProductDiag(const VectorS<M,Tx>& x, const VectorS<M,Td>& d)
  {
    typedef typename promote_Surreal<Tx,Td>::type T;
    T sum = 0.0;
    for (int i = 0; i < M; i++)
        sum += x(i)*d(i)*x(i);

    return sum;
  }

} //namespace DLA
} //namespace tinymat 

#endif //VECTORS_INNERPRODUCT_H
