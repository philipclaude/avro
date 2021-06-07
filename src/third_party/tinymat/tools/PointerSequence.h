// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef POINTERSEQUENCE_H
#define POINTERSEQUENCE_H

#include <algorithm> // rotate

#include "SANSException.h"

namespace tinymat 
{

//----------------------------------------------------------------------------//
// A class that represents a sequence of pointers

template<class T>
class PointerSequence
{
public:
  template< class... Args >
  PointerSequence(const int size, Args&&... args);
  ~PointerSequence();

  // No-copy constructors
  PointerSequence(const PointerSequence&) = delete;
  PointerSequence& operator=(const PointerSequence&) = delete;

  int size() const { return size_; }

        T& operator[](const int i)       { return *data_[i]; }
  const T& operator[](const int i) const { return *data_[i]; }

  // Rotates the last data so it's first. i.e. {d0, d1, d2, d3} -> {d3, d0, d1, d2}
  void rotate()
  {
    // Rotate all the pointers
    std::rotate(data_, data_ +size_-1, data_+size_);
  }

protected:

  T** data_;
  const int size_;
};

template<class T>
template< class... Args >
PointerSequence<T>::PointerSequence(const int size, Args&&... args) : data_(NULL), size_(size)
{
  // There needs to be at least one field
  SANS_ASSERT( size_ >= 0 );

  data_ = new T*[size_];

  for (int i = 0; i < size_; i++ )
    data_[i] = new T(std::forward<Args>(args)...);
}

template<class T>
PointerSequence<T>::~PointerSequence()
{
  for (int i = 0; i < size_; i++ )
    delete data_[i];

  delete [] data_;
}

}

#endif // POINTERSEQUENCE_H
