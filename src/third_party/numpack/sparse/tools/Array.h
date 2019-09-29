// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef ARRAY_H
#define ARRAY_H

#include <iostream>
#include <iomanip>

namespace numpack 
{
namespace SLA
{
//===========================================================================
template<typename T>
class Array
{
public:
  typedef T value_type;
  typedef Array type;

  typedef unsigned const int INT;

public:
  explicit Array(INT n) : size_(n), values_(new T[n]) {}
  Array(const Array&) = delete;
  ~Array() { delete [] values_; }

  inline T& operator[](INT i) const { return values_[i]; }

  Array<T>& operator=(T const &val)
  {
    for (unsigned int i = 0; i < size_; ++i)
      values_[i] = val;

    return *this;
  }

  template<typename T_>
  friend std::ostream& operator<<(std::ostream& out, Array<T_> const &A); // output

private:
  unsigned const int size_; //size of the array
  T *values_;               //Storage space for the array

};

//=============================================================================
template<typename T>
std::ostream& operator<<(std::ostream& out, Array<T> const &A) // output
{
  unsigned int i;
  for (i = 0; i < A.size_; i++)
    out << std::setw(15) << std::setprecision(4) << A[i];

  return out;
}

} // namespace SLA
} // namespace numpack 

#endif // ARRAY_H
