// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef PETSC_VECTORSIZE_IMPL_HACK_H
#define PETSC_VECTORSIZE_IMPL_HACK_H

#include "tools/SANSnumerics.h"

namespace numpack 
{
namespace DLA
{
template<int M, class T>
class VectorS;

template<class T>
class VectorDView;
}

namespace SLA
{

template<class T>
class SparseVector;


//TODO: This is a complete HACK to get the continuous global mapping
template< class Vector_type >
struct StateVectorSize;

template<>
struct StateVectorSize< SLA::SparseVector<Real> > { static const int M = 1; };

template<int M_>
struct StateVectorSize< SLA::SparseVector<DLA::VectorS<M_,Real>> > { static const int M = M_; };

template<>
struct StateVectorSize< DLA::VectorDView<SLA::SparseVector<Real>> > { static const int M = 1; };

template<int M_>
struct StateVectorSize< DLA::VectorDView<SLA::SparseVector<DLA::VectorS<M_,Real>>> > { static const int M = M_; };


} //namespace SLA
} //namespace numpack 

#endif //PETSC_VECTORSIZE_IMPL_HACK_H
