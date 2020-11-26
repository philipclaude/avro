// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXS_EIGEN_H
#define MATRIXS_EIGEN_H

#include "MatrixS_Type.h"

//Eigen value/vector interface for dense.

namespace tinymat 
{
namespace DLA
{
  template<int M, class T>
  struct EigenSystemPair;

  template< int M, class T >
  inline void
  EigenValues( MatrixS<M,M,T>& A, VectorS<M,T>& wr, VectorS<M,T>& wi );

  template< int M, class T >
  inline void
  EigenVectors( MatrixS<M,M,T>& A, MatrixS<M,M,T>& vl, MatrixS<M,M,T>& vr );

  template< int M, class T >
  inline void
  EigenSystem( MatrixS<M,M,T>& A, VectorS<M,T>& wr, VectorS<M,T>& wi, MatrixS<M,M,T>& vl, MatrixS<M,M,T>& vr );


  template< int M, class T >
  void
  EigenValues( const MatrixSymS<M,T>& A, VectorS<M,T>& L );

  template< int M, class T >
  void
  EigenVectors( const MatrixSymS<M,T>& A, MatrixS<M,M,T>& E );

  template< int M, class T >
  void
  EigenSystem( const MatrixSymS<M,T>& A, VectorS<M,T>& L, MatrixS<M,M,T>& E );


  template<int M, class T>
  struct EigenSystemPair
  {
    // cppcheck-suppress noExplicitConstructor
    EigenSystemPair(const MatrixSymS<M,T>& A) { EigenSystem(A, L, E); }
    EigenSystemPair(const EigenSystemPair<M,T>& LE ) : L(LE.L), E(LE.E) {}
    EigenSystemPair& operator=(const EigenSystemPair<M,T>& LE ) { L = LE.L; E = LE.E; return *this; }
    VectorS<M,T> L;
    MatrixS<M,M,T> E;
  };

  template< int M, class T >
  void
  EigenSystem( const MatrixSymS<M,T>& A, EigenSystemPair<M,T>& LE ) { EigenSystem(A, LE.L, LE.E); }

} //namespace tinymat 
} //namespace DLA


#endif //MATRIXS_EIGEN_H
