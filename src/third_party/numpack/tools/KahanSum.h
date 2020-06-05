// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef KAHANSUM_H
#define KAHANSUM_H

// Kahan summation reduces numerical error when summing a sequence
// of finite precision floating point numbers
//
// http://codereview.stackexchange.com/questions/56532/kahan-summation
//
namespace numpack 
{

template<class T>
class KahanSum
{
public:

  KahanSum() : sum_(0), running_error_(0) {}
  // cppcheck-suppress noExplicitConstructor
  KahanSum(const T& val) : sum_(val), running_error_(0) {}

  T operator+=(const T& val)
  {
#ifdef __INTEL_COMPILER
    // The intel compiler tends to remove the presidence when computing the rounding error...
    T difference = val - running_error_;
    T temp = sum_ + difference;
    running_error_ = diff(temp - sum_) - difference;
    sum_ = temp;
#else
    T difference = val - running_error_;
    T temp = sum_ + difference;
    running_error_ = (temp - sum_) - difference;
    sum_ = temp;
#endif

    return sum_;
  }

  T operator=(const T val)
  {
    sum_ = val;
    running_error_ = 0;
    return sum_;
  }

  T operator/=(const T val)
  {
    sum_ = sum_/val;
    return sum_;
  }

  operator const T() const { return sum_; }

protected:
#ifdef __INTEL_COMPILER
  __attribute__((noinline))
  T diff( const T val ) const;
#endif

  T sum_;
  T running_error_;
};

#ifdef __INTEL_COMPILER
// Intel creates a warning if the method is defined inside the body
template<class T>
__attribute__((noinline))
T KahanSum<T>::diff( const T val ) const
{
  return val;
}
#endif

}

#endif // KAHANSUM_H
