// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXD_TYPE_H
#define MATRIXD_TYPE_H

namespace numpack
{
namespace DLA
{

//Base class to distinguish what types belong to the dense matrix family
//The bool flag indicates if recursive function expressions should be used
template< class Derived, bool useRF >
class MatrixDType
{
public:
  //A convenient method for casting to the derived type
  inline const Derived& cast() const { return static_cast<const Derived&>(*this); }

  int m() const { return cast().m(); }
  int n() const { return cast().n(); }
  int size() const { return m()*n(); }
};

//Forward declarations
template< class T >
class MatrixDView;

template< class T >
class MatrixD;

template <class T>
class MatrixSymD;

template< class MatrixL >
class MatrixDTuple;

template< class T >
class VectorDView;

template< class T >
class VectorD;

template<class ExprL, class ExprR>
class OpMulD;

class DenseMatrixSize;
class DenseVectorSize;

//Represents the transposed view of a matrix
template<class T>
class MatrixDTranspose;

// Dummy non-zero pattern for Dense matrix for symmetry with SparseLinAlg
template< class T >
class DenseNonZeroPattern;

}
}

#endif //MATRIX_TYPE_H
