// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef sparse_TYPE_H
#define sparse_TYPE_H

namespace numpack 
{
namespace SLA
{

//Base class to distinguish what types belong to the sparse linear algebra family
//useRF is a flag to indicate the expression includes a matmul
template< class Derived, bool useRF >
class sparseType
{
public:
  //A convenient method for casting to the derived type
  inline const Derived& cast() const { return static_cast<const Derived&>(*this); }

  int m() const { return cast().m(); }
};


// Forward declarations
template<class TM>
class SparseMatrix_CRS;

template<class TV>
class SparseVector;

class SparseMatrixSize;

//Represents the transpose of a matrix
template<class TM_>
class SparseMatrix_CRS_Transpose;

//Non-zero pattern for a sparse matrix
template< class TM >
class SparseNonZeroPattern;

//Provides a transposed view of a sparse non-zero pattern
template< class TM >
class SparseNonZeroPattern_Transpose;
}
}

#endif //sparse_TYPE_H
