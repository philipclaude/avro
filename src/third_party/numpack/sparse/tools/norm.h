// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef sparse_NORM_H
#define sparse_NORM_H

#include "tools/SANSnumerics.h"

#include "numpack/sparse/SparseVector.h"
#include "numpack/sparse/tools/sparse_Scalar.h"

#include "numpack/dense/tools/norm.h"

namespace numpack 
{

inline Real sum_pow(const Real& a, const unsigned int& p);

template<int N>
inline SurrealS<N> sum_pow(const SurrealS<N>& a, const unsigned int& p);


template<class T>
inline typename Scalar<T>::type sum_pow(const SLA::SparseVector<T>& a, const unsigned int& p)
{
  const unsigned int size = a.m();

  typename Scalar<T>::type val = 0;
  for (unsigned int i = 0; i < size; i++)
    val += sum_pow(a[i],p);

  return val;
}

template<class T>
inline typename Scalar<T>::type norm(const SLA::SparseVector<T>& a, const unsigned int& p)
{
  const unsigned int size = a.m();

  typename Scalar<T>::type val = 0;
  for (unsigned int i = 0; i < size; i++)
    val += sum_pow(a[i],p);

  return pow(val,1./Real(p));
}

}

#endif //sparse_NORM_H
