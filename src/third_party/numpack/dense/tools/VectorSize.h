// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef DLA_VECTORSIZE_H
#define DLA_VECTORSIZE_H

#include "tools/SANSnumerics.h"
#include "numpack/dense/static/MatrixS_Type.h"

// forward declaration
template<int N, class T>
class SurrealS;

namespace numpack 
{
namespace DLA
{

template<class T>
struct VectorSize;

// Gives the size of a VectorS or Real

template<int M_, class T> struct VectorSize< VectorS<M_,T> > { static const int M = M_; };

template<> struct VectorSize<Real> { static const int M = 1; };
template<> struct VectorSize<int > { static const int M = 1; };
template<int N, class T> struct VectorSize<SurrealS<N,T>> { static const int M = 1; };

}
}


#endif //DLA_VECTORSIZE_H
