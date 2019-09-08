// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef POINTERARRAY_H
#define POINTERARRAY_H

#include <algorithm> // rotate

#include "SANSException.h"

namespace SANS
{

//----------------------------------------------------------------------------//
// A class that represents an array of pointers

template<class T>
class PointerArray
{
public:
  PointerArray() : size_(0), data_(nullptr) {}
  explicit PointerArray(const int size);
  ~PointerArray();

  // No-copy constructors
  PointerArray(const PointerArray&) = delete;
  PointerArray& operator=(const PointerArray&) = delete;

  int size() const { return size_; }

        T*& operator[](const int i)       { return data_[i]; }
  const T*  operator[](const int i) const { return data_[i]; }

  // Rotates the last data so it's first. i.e. {d0, d1, d2, d3} -> {d3, d0, d1, d2}
  void rotate()
  {
    // Rotate all the pointers
    std::rotate(data_, data_ +size_-1, data_+size_);
  }

  // Resizes the array of pointers keeping old data that still fits
  void resize(const int size);

protected:

  int size_;
  T** data_;
};

template<class T>
PointerArray<T>::PointerArray(const int size) : size_(size), data_(NULL)
{
  // There needs to be at least one field
  SANS_ASSERT( size_ >= 0 );

  if ( size_ == 0 ) return;

  data_ = new T*[size_];

  for (int i = 0; i < size_; i++ )
    data_[i] = nullptr;
}

template<class T>
PointerArray<T>::~PointerArray()
{
  for (int i = 0; i < size_; i++ )
    delete data_[i];

  delete [] data_;
}

template<class T>
void
PointerArray<T>::resize(const int size)
{
  SANS_ASSERT( size >= 0 );

  // delete all pointers for a size of 0
  if ( size == 0 )
  {
    for (int i = 0; i < size_; i++)
      delete data_[i];

    delete [] data_;
    data_ = nullptr;
    size_ = 0;

    return;
  }

  T** olddata = data_;
  int oldsize = size_;
  size_ = size;

  // delete only those groups larger than the new size
  for (int i = size; i < oldsize; i++)
    delete data_[i];

  // allocate the new memory and initialize to null
  data_ = new T*[size];

  for (int i = 0; i < size_; i++)
    data_[i] = nullptr;

  // copy over the old data pointers
  for (int i = 0; i < std::min(oldsize, size_); i++)
    data_[i] = olddata[i];

  // remove the old array of pointers
  delete [] olddata;
}

}

#endif // POINTERSEQUENCE_H
