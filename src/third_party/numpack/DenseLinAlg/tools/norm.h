// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef DENSELINALG_NORM_H
#define DENSELINALG_NORM_H

#include <cmath> // pow

#include "tools/SANSnumerics.h"

#include "numpack/DenseLinAlg/StaticSize/MatrixS_Type.h"
#include "numpack/DenseLinAlg/DynamicSize/MatrixD_Type.h"

#include "tools/SANSTraitsScalar.h"

#include "numpack/types/SurrealS_Type.h"

namespace numpack 
{

namespace SLA
{
template<class TV>
class SparseVector;
}

template<class T>
inline typename Scalar<T>::type sum_pow(const SLA::SparseVector<T>& a, const unsigned int& p);

template<class T>
inline typename Scalar<T>::type norm(const SLA::SparseVector<T>& a, const unsigned int& p);


inline Real sum_pow(const Real& a, const unsigned int& p) { return pow(a,p); }

template<int N>
inline SurrealS<N> sum_pow(const SurrealS<N>& a, const unsigned int& p) { return pow(a,p); }

template<int M, class T>
inline typename Scalar<T>::type sum_pow(const DLA::MatrixS<M,1,T>& a, const unsigned int& p)
{
  typename Scalar<T>::type val = 0;
  for (int i = 0; i < M; i++)
    val += sum_pow(a(i,0),p);

  return val;
}

template<int M, class T>
inline typename Scalar<T>::type norm(const DLA::MatrixS<M,1,T>& a, const unsigned int& p)
{
  typename Scalar<T>::type val = 0;
  for (int i = 0; i < M; i++)
    val += sum_pow(a(i,0),p);

  return pow(val,1./p);
}

template<class T>
inline typename Scalar<T>::type sum_pow(const DLA::VectorDView<T>& a, const unsigned int& p)
{
  const int size = a.size();

  typename Scalar<T>::type val = 0;
  for (int i = 0; i < size; i++)
    val += sum_pow(a[i],p);

  return val;
}

template<class T>
inline typename Scalar<T>::type norm(const DLA::VectorDView<T>& a, const unsigned int& p)
{
  const int size = a.size();

  typename Scalar<T>::type val = 0;
  for (int i = 0; i < size; i++)
    val += sum_pow(a[i],p);

  return pow(val,1./p);
}

}

#endif //DENSELINALG_DOT_H
