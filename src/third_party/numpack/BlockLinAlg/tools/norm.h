// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef BLOCKLINALG_NORM_H
#define BLOCKLINALG_NORM_H

#include "tools/SANSnumerics.h"

#include "numpack/SparseLinAlg/SparseVector.h"
#include "numpack/SparseLinAlg/tools/SparseLinAlg_Scalar.h"

#include "numpack/DenseLinAlg/tools/norm.h"

#include "numpack/types/PromoteSurreal.h"

namespace numpack 
{

template<class T>
inline typename Scalar<T>::type norm(BLA::VectorBlock_2<DLA::VectorD<SLA::SparseVector<T> >,
                                                        DLA::VectorD<SLA::SparseVector<T> > >& a,
                                     const unsigned int& p)
{

  return pow(pow(norm(a.v0, p), p) + pow(norm(a.v1, p), p),1./Real(p));
}

template<class T0, class T1>
inline typename promote_Surreal<typename Scalar<T0>::type, typename Scalar<T1>::type>::type
norm(BLA::VectorBlock_2<DLA::VectorD<SLA::SparseVector<T0> >,
                        DLA::VectorD<SLA::SparseVector<T1> > >& a,
     const unsigned int& p)
{
//  SANS_DEVELOPER_EXCEPTION("This is not implemented");
  return pow(pow(norm(a.v0, p), p) + pow(norm(a.v1, p), p),1./Real(p));
}


template<class T0, class T1, class T2, class T3>
inline typename promote_Surreal<typename Scalar<T0>::type,
                                typename Scalar<T1>::type,
                                typename Scalar<T2>::type,
                                typename Scalar<T3>::type>::type
norm(BLA::VectorBlock_4<DLA::VectorD<SLA::SparseVector<T0> >,
                        DLA::VectorD<SLA::SparseVector<T1> >,
                        DLA::VectorD<SLA::SparseVector<T2> >,
                        DLA::VectorD<SLA::SparseVector<T3> > >& a,
     const unsigned int& p)
{
//  SANS_DEVELOPER_EXCEPTION("This is not implemented");
  return pow( pow(norm(a.v0, p), p) + pow(norm(a.v1, p), p) + pow(norm(a.v2, p), p) + pow(norm(a.v3, p), p), 1./Real(p));
}

}

#endif //BLOCKLINALG_NORM_H
