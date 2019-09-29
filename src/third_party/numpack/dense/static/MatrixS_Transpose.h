// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXS_TRANSPOSE_H
#define MATRIXS_TRANSPOSE_H

#include <stdint.h> // uintptr_t

#include "tools/SANSnumerics.h"

#include "MatrixS_Type.h"
#include "numpack/Transpose.h"

namespace numpack 
{
namespace DLA
{

//Represents the transpose of a matrix
template<int M_, int N_, class T>
class MatrixSTranspose : public MatrixSType< MatrixSTranspose<M_,N_,T>, true, true >
{

public:
//  typedef T Ttype;
  typedef typename TransposeTraits<T>::type Ttype;
  typedef typename TransposeViewTraits<T>::type Ttransview;

  //Transpose
  static const int M = N_;
  static const int N = M_;

  explicit MatrixSTranspose( MatrixS<M_,N_,T>& A ) : A(A) {}

  //Operator to access the matrix values
  inline       Ttransview operator()(const int i, const int j)       { return Transpose(A(j,i)); }
  inline const Ttransview operator()(const int i, const int j) const { return Transpose(A(j,i)); }

  // Lazy expression operations
  template<class Scalar>
  inline void value(const Scalar& sgn, MatrixS<M,N,Ttype>& res) const
  {
    for (int i = 0; i < M; ++i)
      for (int j = 0; j < N; ++j)
        res(i,j) = sgn*(*this)(i,j);
  }
  template<class Scalar>
  inline void plus(const Scalar& sgn, MatrixS<M,N,Ttype>& res) const
  {
    for (int i = 0; i < M; ++i)
      for (int j = 0; j < N; ++j)
        res(i,j) += sgn*(*this)(i,j);
  }

  //A unique ID for the matrix so the matrix can be uniquely identified
  uintptr_t ID() const { return A.ID(); }

private:
  MatrixS<M_,N_,T>& A; //The matrix to transpose
};

} //namespace DLA
} //namespace numpack 



#endif //MATRIXS_TRANSPOSE_H
