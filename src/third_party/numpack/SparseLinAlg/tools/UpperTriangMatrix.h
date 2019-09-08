// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef UPPER_TRIANG_MATRIX_H
#define UPPER_TRIANG_MATRIX_H

#include "tools/SANSnumerics.h"

#include <iostream>

namespace numpack 
{
namespace SLA
{
//===========================================================================
// A class for efficiently storing an upper triangular matrix
class UpperTriangMatrix
{
public:

  typedef unsigned const int INT;
  explicit UpperTriangMatrix(INT n) : n_(n), size_(n*(n+1)/2), values_(new Real[size_]) { *this = 0;}
  UpperTriangMatrix( const UpperTriangMatrix& M ) : n_(M.n_), size_(M.size_), values_(new Real[size_])
  {
    for (unsigned int i = 0; i < size_; ++i)
      values_[i] = M.values_[i];
  }
  UpperTriangMatrix& operator=( const UpperTriangMatrix& ) = delete;

  ~UpperTriangMatrix() { delete [] values_; }

  inline Real& operator()(INT i, INT j) const { return values_[size_ - (n_-i)*((n_-i)+1)/2 - i + j]; }

  inline UpperTriangMatrix& operator=(Real const &val)
  {
    for (unsigned int i = 0; i < size_; ++i)
      values_[i] = val;

    return *this;
  }

  friend std::ostream& operator<<(std::ostream& out, UpperTriangMatrix const &TriM);

private:
  unsigned const int n_; //Number of rows/columns
  unsigned const int size_; //Total size of the memory
  Real *values_; //Storage space for the triangular matrix
};

//=============================================================================
std::ostream& operator<<(std::ostream& out, UpperTriangMatrix const &TriM); // output


} // namespace SLA
} // namespace numpack 

#endif //UPPER_TRIANG_MATRIX_H
