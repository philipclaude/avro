// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef VECTORD_CLASS_H
#define VECTORD_CLASS_H

// Vector class with run-time size

#include <iostream>

#include "tools/SANSTraitsPOD.h"
#include "tools/SANSException.h"

#include "MatrixD.h"

namespace numpack 
{
namespace DLA
{


//----------------------------------------------------------------------------//
// MatrixDView:  Linear Algebra Dense Matrix
//
// Operators with Lazy Recursive Functions
//
// MatrixDView does not allocate/deallocate the memory
//----------------------------------------------------------------------------//

template< class T >
class VectorDView : public MatrixDView<T>
{
public:
  typedef T node_type;

protected:
  // No default constructor
  VectorDView() : MatrixDView<T>(NULL, 0, 1) {}

public:
  explicit VectorDView( T v0[], const int m ) : MatrixDView<T>(v0, m, 1) {}
  explicit VectorDView( T v0[], const int m, const int stride ) : MatrixDView<T>(v0, m, 1, stride) {}
  VectorDView( T v0[], const int m, const typename POD<T>::type& v ) : MatrixDView<T>(v0, m, 1, v) {}
  VectorDView( const VectorDView& z ) : MatrixDView<T>(z) {}
  ~VectorDView() {}

  //Extract a sub-vector
  VectorDView sub( const int i, const int m ) const
  {
    return VectorDView(v_ + i*stride_, m, stride_ );
  }
  using MatrixDView<T>::sub;

  // access operators
        T& operator[]( int n )       { return v_[n*stride_]; }
  const T& operator[]( int n ) const { return v_[n*stride_]; }
        T& operator()( int n )       { return v_[n*stride_]; }
  const T& operator()( int n ) const { return v_[n*stride_]; }

  // assignment
  VectorDView& operator=( const VectorDView& v )           { MatrixDView<T>::operator=(v); return *this; }
  VectorDView& operator=( const T& v )                     { MatrixDView<T>::operator=(v); return *this; }
  VectorDView& operator=( const typename POD<T>::type& v ) { MatrixDView<T>::operator=(v); return *this; }

  // Assignment operator with expression templates
  template<class Expr, bool useRF>
  VectorDView& operator= ( const MatrixDType<Expr, useRF>& r ) { MatrixDView<T>::operator=(r); return *this; }

protected:
  using MatrixDView<T>::v_;      // value
  using MatrixDView<T>::m_;      // m dimension of the vector
  using MatrixDView<T>::stride_; // Stride for when the memory is larger than m
};


//----------------------------------------------------------------------------//
// VectorD<T>:  vector w/ data elements of type T
//
// VectorD allocates/deallocates memory
//
// template parameters:
//   T            array data type (e.g. double)
//----------------------------------------------------------------------------//

template <class T>
class VectorD : public VectorDView<T>
{
public:
  typedef T node_type;
  typedef typename MatrixDView<T>::size_type size_type;

  friend class VectorD< DenseVectorSize >;
  friend class VectorD< VectorD<T> >;
protected:
  // Needed create MatrixD< MatrixD<T> >
  VectorD() : VectorDView<T>() {}

  //Used in DLA::VectorD to allocate an array of VectorD
  VectorD& operator=(const DenseVectorSize& s)
  {
    delete [] v_; v_ = NULL;
    v_ = new T[s.m()];
    m_ = s.m();
    return *this;
  }
public:

  explicit VectorD( const int m ) : VectorDView<T>( new T[m], m ) {}
  VectorD( const VectorD& v   ) : VectorDView<T>(new T[v.m_], v.m_) { init_array_copy(v.m_, v.v_, v_); MatrixDView<T>::operator=(v); }
  VectorD( const size_type& s ) : VectorDView<T>(new T[s.m_], s.m_) { init_array_copy(s.m_, s.v_, v_); }
  VectorD( const DenseVectorSize& s ) : VectorDView<T>(new T[s.m()], s.m()) {}
  VectorD( const int m, const T& v0 ) : VectorDView<T>(new T[m], m) { MatrixDView<T>::operator=(v0); }
  VectorD( const int m, const typename numpack::POD<T>::type& v0 ) : VectorDView<T>(new T[m], m, v0) { MatrixDView<T>::operator=(v0); }
  VectorD( const int m, const T v0[] ) : VectorDView<T>(new T[m], m)
  {
    for (int i = 0; i < m_; i++ ) v_[i] = v0[i];
  }
  // Move constructor just takes the memory
  VectorD( VectorD&& v ) : VectorDView<T>(v.v_, v.m_) { v.v_ = NULL; }

  template<class U>
  VectorD( const VectorD<U>& v );
  template<class U>
  VectorD( const std::initializer_list<U>& s );
#ifdef __INTEL_COMPILER
  VectorD( const std::initializer_list< std::initializer_list<int> >& s );
#else
  template<class U>
  VectorD( const std::initializer_list< std::initializer_list<U> >& s );
#endif
  ~VectorD() { delete [] v_; }

  // access operators
  // duplicated here to remove arithemtic operation on access
        T& operator[]( int n )       { return v_[n]; }
  const T& operator[]( int n ) const { return v_[n]; }
        T& operator()( int n )       { return v_[n]; }
  const T& operator()( int n ) const { return v_[n]; }

  //This allows a vector to be assigned an expression as it is constructed, i.e.
  //VectorD<T> C = A + B;
  template<class Expr, bool useRF>
  VectorD( const MatrixDType<Expr, useRF>& r ) : VectorDView<T>( new T[r.cast().m()], r.cast().m() ) { MatrixDView<T>::operator=(r); }

  void resize( const size_type& s )
  {
    delete [] v_; v_ = NULL;
    v_ = new T[s.m_];
    m_ = s.m_;
    init_array_copy(s.m_, s.v_, v_);
  }

  void resize( const int m )
  {
    delete [] v_; v_ = NULL;
    v_ = new T[m];
    m_ = m;
  }

  // assignment
  VectorD& operator=( const VectorD& v )                                                  { MatrixDView<T>::operator=(v); return *this; }
  VectorD& operator=( const T& s )                                                        { MatrixDView<T>::operator=(s); return *this; }
  VectorD& operator=( const typename numpack::POD<T>::type& s )                              { MatrixDView<T>::operator=(s); return *this; }
  VectorD& operator=( const std::initializer_list<T>& s );

  template<class Expr, bool useRF>
  VectorD& operator=( const MatrixDType<Expr, useRF>& r)
  {
    // Allow the vector to be resized if it's NULL
    if (v_ == NULL) resize(r.cast().m());

    MatrixDView<T>::operator=(r);
    return *this;
  }

  VectorD& operator=( VectorD&& v )
  {
    if (v_ == NULL)
    {
      // Simply steal the memory
      m_ = v.m();
    }
    else
    {
      // Check the size match before stealing the memory
      SANS_ASSERT( m_ == v.m_ );
      delete [] v_; v_ = NULL;
    }

    // Steal the memory
    v_ = v.v_;
    v.v_ = NULL;

    return *this;
  }

  void dump( int indentSize=0, std::ostream& out = std::cout ) const;

protected:
  using VectorDView<T>::v_;
  using VectorDView<T>::m_;
};

template <class T>
template <class U>
VectorD<T>::VectorD( const VectorD<U>& v ) : VectorDView<T>(new T[v.m()], v.m())
{
  init_array_copy(m_, &v[0], v_);
  for (int i = 0; i < v.m(); ++i)
    v_[i] = v[i];
}

template <class T>
template <class U>
VectorD<T>::VectorD( const std::initializer_list<U>& s ) : VectorDView<T>(new T[s.size()], s.size())
{
  init_array_copy(m_, s.begin(), v_);
  int n = 0;
  for (auto i = s.begin(); i != s.end(); ++i)
    v_[n++] = *i;
}

#if __INTEL_COMPILER
//Maybe someday the intel compiler will get fixed and we won't need any of this mess...
template <class T>
VectorD<T>::VectorD( const std::initializer_list< std::initializer_list<int> >& s ) : VectorDView<T>(new T[s.size()], s.size())
{
  int n = 0;
  for (auto i = s.begin(); i != s.end(); ++i)
  {
    initializer_list_assign<T>(v_[n],*i);
    n++;
  }
}
#else
template <class T>
template <class U>
VectorD<T>::VectorD( const std::initializer_list< std::initializer_list<U> >& s ) : VectorDView<T>(new T[s.size()], s.size())
{
  init_array_copy(m_, s.begin(), v_);
  int n = 0;
  for (auto i = s.begin(); i != s.end(); ++i)
    v_[n++] = *i;
}
#endif

template <class T>
inline VectorD<T>&
VectorD<T>::operator=( const std::initializer_list<T>& s )
{
  SANS_ASSERT(s.size() == static_cast<std::size_t>(this->size()));
  int n = 0;
  for (typename std::initializer_list<T>::iterator i = s.begin(); i != s.end(); ++i)
    v_[n++] = *i;
  return *this;
}

// debug dump of private data
template <class T>
void
VectorD<T>::dump( int indentSize, std::ostream& out ) const
{
  std::string indent(indentSize, ' ');
  out << indent << "VectorD<T>:" << std::endl;
#if 1
  out << indent << "  v_ = ";
  for (int n = 0; n < m_; n++)
    out << v_[n] << " ";
  out << std::endl;
#else     // only works for class T with member function dump()
  for (int n = 0; n < N; n++)
  {
    out << indent << "  v_[" << n << "] = ";
    data[n].dump(indentSize, out);
  }
#endif
}

// I/O

template <class T>
std::ostream&
operator<<( std::ostream& out, const VectorD<T>& v )
{
  for (int i = 0; i < v.m()-1; i++)
    out << v[i] << ", " << std::endl << std::endl;
  if (v.m() > 0 )
    out << v[v.m()-1];
  return out;
}

//These member functions of MatrixDView explicitly depend on VectorDView, hence they must be defined here rather than in MatrixD.h
template<class T>
VectorDView<T>
MatrixDView<T>::subcol( const int i, const int j, const int n ) const { return VectorDView<T>(v_ + i*stride_ + j, n, stride_ ); }

template<class T>
VectorDView<T>
MatrixDView<T>::subcol( const int i, const int j, const int n )       { return VectorDView<T>(v_ + i*stride_ + j, n, stride_ ); }



//----------------------------------------------------------------------------//
// VectorD<DenseVectorSize>:  vector w/ data elements of typeDenseVectorSizeT
//
// VectorD allocates/deallocates memory
//----------------------------------------------------------------------------//

template <>
class VectorD< DenseVectorSize > : public VectorDView< DenseVectorSize >
{
public:
  typedef DenseVectorSize node_type;
  typedef DenseVectorSize T;
  typedef typename MatrixDView<DenseVectorSize>::size_type size_type;

public:

  explicit VectorD( const int m ) : VectorDView<T>( new T[m], m ) {}

  // Copy constructor
  VectorD( const VectorD& v ) : VectorDView<T>(new T[v.m_], v.m_) { MatrixDView<T>::operator=(v); }

  // Move constructor just takes the memory
  VectorD( VectorD&& v ) : VectorDView<T>(v.v_, v.m_) { v.v_ = NULL; }

  VectorD( const std::initializer_list< int >& s );

  ~VectorD() { delete [] v_; }

  // access operators
  // duplicated here to remove arithemtic operation on access
        T& operator[]( int n )       { return v_[n]; }
  const T& operator[]( int n ) const { return v_[n]; }
        T& operator()( int n )       { return v_[n]; }
  const T& operator()( int n ) const { return v_[n]; }

  // assignment
  VectorD& operator=( const VectorD& v )                     { MatrixDView<T>::operator=(v); return *this; }
  VectorD& operator=( const T& s )                           { MatrixDView<T>::operator=(s); return *this; }

  void dump( int indentSize=0, std::ostream& out = std::cout ) const;

protected:
  using VectorDView<T>::v_;
  using VectorDView<T>::m_;
};

inline VectorD<DenseVectorSize>::VectorD( const std::initializer_list< int >& s )
  : VectorDView<DenseVectorSize>(new DenseVectorSize[s.size()], s.size())
{
  int n = 0;
  for (auto i = s.begin(); i != s.end(); ++i)
    v_[n++] = *i;
}

} //namespace DLA
} //namespace numpack 

#endif // VECTORD_CLASS_H
