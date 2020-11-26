// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef INIT_ARRAY_COPY_H
#define INIT_ARRAY_COPY_H

namespace tinymat 
{
  //An empty function that can be used to call special initialization
  //methods of when creating copies of classes
  template<class T0, class T1>
  inline void init_array_copy( const int size, const T0* orig, T1* copy ) {}

  namespace DLA
  {
    template<class TM>
    class MatrixD;

    template<class TV>
    class VectorD;
  }

  // Specialization for dense matrix
  template<class TM>
  inline
  void init_array_copy( const int size, const DLA::MatrixD<TM>* orig, DLA::MatrixD<TM>* copy )
  {
    for ( int i = 0; i < size; i++ )
      copy[i].resize(orig[i].size());
  }

  // Specialization for dense vectors
  template<class TV>
  inline
  void init_array_copy( const int size, const DLA::VectorD<TV>* orig, DLA::VectorD<TV>* copy )
  {
    for ( int i = 0; i < size; i++ )
      copy[i].resize(orig[i].size());
  }

  namespace SLA
  {
    template<class TV>
    class SparseVector;
  }

  //Specialization for copies of SparseVector
  template<class TV>
  inline
  void init_array_copy( const int size, const SLA::SparseVector<TV>* orig, SLA::SparseVector<TV>* copy )
  {
    for ( int i = 0; i < size; i++ )
      copy[i].resize(orig[i].size());
  }
}

#endif //INIT_ARRAY_COPY_H
