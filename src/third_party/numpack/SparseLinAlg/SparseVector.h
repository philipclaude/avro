// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SPARSEVECTOR_SCALAR_H
#define SPARSEVECTOR_SCALAR_H

#include <cmath> // sqrt

#include "tools/SANSnumerics.h"
#include "tools/SANSException.h"
#include "tools/SANSTraitsInitListAssign.h"

#include "SparseLinAlg_Type.h"
#include "numpack/DenseLinAlg/DynamicSize/MatrixD.h"
#include "numpack/DenseLinAlg/DynamicSize/VectorD.h"
#include "numpack/DenseLinAlg/tools/dot.h"

namespace numpack 
{

#ifdef __INTEL_COMPILER
namespace SLA
{
//Forward declare
template<class TV>
class SparseVector;
}
//
// This is completely unnessary if the intel compiler could use templated initializer_list functions....
//
template<class TV>
struct initializer_list_assign< SLA::SparseVector< TV > >
{
  initializer_list_assign(SLA::SparseVector< TV >& val, const std::initializer_list<int>& s ) { val = s; }
};
#endif

namespace SLA
{

struct SparseVectorSize
{
  SparseVectorSize() : m_(-1) {}
  explicit SparseVectorSize(const int m) : m_(m) {}
  SparseVectorSize(const SparseVectorSize& size) : m_(size.m_) {}
  SparseVectorSize& operator=(const int m ) { m_ = m; return *this; }
  SparseVectorSize& operator=(const SparseVectorSize& s ) { m_ = s.m_; return *this; }
  //operator int() const { return m_; }
  int m() const {return m_;}
  int m_;
};


template<class TV>
class SparseVector : public SparseLinAlgType< SparseVector<TV>, false >
{
public:
  typedef TV Ttype;

  typedef SparseVectorSize size_type;

  explicit SparseVector( const unsigned int m ) : m_(m), values_( new TV[m_] ) {}
  // cppcheck-suppress noExplicitConstructor
  SparseVector( const size_type& size ) : m_(size.m_), values_( new TV[m_] ) {}
  SparseVector( const SparseVector& v ) : m_(v.m_), values_( new TV[m_] ) { *this = v; }

  // move constructor takes the memory
  SparseVector( SparseVector&& v ) : m_(v.m_), values_( v.values_ ) { v.values_ = NULL; }

  template<class Expr, bool useRF> // cppcheck-suppress noExplicitConstructor
  SparseVector( const SparseLinAlgType<Expr, useRF>& r ) : m_(r.m()), values_( new TV[m_] ) { *this = r; }


  friend class DLA::VectorD< SparseVector<TV> >;
  friend class DLA::VectorD< SparseVector<size_type> >;
#ifdef __INTEL_COMPILER
  friend class initializer_list_assign< SparseVector<TV> >;
#endif
protected:
  SparseVector() : m_(0), values_( 0 ) {}
  SparseVector& operator=( const std::initializer_list<int>& s )
  {
    SANS_ASSERT(s.size() == 1);
    resize(*s.begin());
    return *this;
  }

  SparseVector& operator=( const size_type& size )
  {
    resize(size.m_);
    return *this;
  }
public:

  ~SparseVector() { delete [] values_; }

  //Changes the size of the vector by removing all the previous data
  void resize( const unsigned int m )  { delete [] values_; m_ = m;       values_ = new TV[m_]; }
  void resize( const size_type& size ) { delete [] values_; m_ = size.m_; values_ = new TV[m_]; }

  int m() const { return m_; }
  int n() const { return 1; }
  size_type size() const { return size_type(m_); }

  SparseVector& operator=(const Real v);
  SparseVector& operator=(const SparseVector& v);
  SparseVector& operator/=(const Real v);
  SparseVector& operator*=(const Real v);

  //L2 norm
  Real norm() const
  {
    Real norm = 0;
    for (int i = 0; i < m_; i++)
      norm += numpack::dot( (*this)[i],(*this)[i] );

    return sqrt(norm);
  }

  //Dot product with another vector
  Real dot( const SparseVector& v ) const
  {
    SANS_ASSERT( m_ == v.m() );

    Real dot = 0;
    for ( int i = 0; i < m_; i++)
      dot += numpack::dot( (*this)[i], v[i] );

    return dot;
  }

  //Access operators
        TV& operator[](const unsigned int i)       { return values_[i]; }
  const TV& operator[](const unsigned int i) const { return values_[i]; }
        TV& operator()(const unsigned int i)       { return values_[i]; }
  const TV& operator()(const unsigned int i) const { return values_[i]; }

  //Lazy evaluation functions
  void value(const Real sgn, SparseVector& b) const;
  void plus(const Real sgn, SparseVector& b) const;

  //Lazy assignment operators
  template< class Expr > inline SparseVector& operator=( const SparseLinAlgType<Expr, true>& r );
  template< class Expr > inline SparseVector& operator+=( const SparseLinAlgType<Expr, true>& r );
  template< class Expr > inline SparseVector& operator-=( const SparseLinAlgType<Expr, true>& r );

  template< class Expr > inline SparseVector& operator=( const SparseLinAlgType<Expr, false>& r );
  template< class Expr > inline SparseVector& operator+=( const SparseLinAlgType<Expr, false>& r );
  template< class Expr > inline SparseVector& operator-=( const SparseLinAlgType<Expr, false>& r );

protected:
  int m_;      //Number of rows in the vector
  TV *values_; //Values of the vector
};

} //namespace SLA
} //namespace numpack 


#define TV_SPEC TV
#include "SparseVector_impl.h"
#undef TV_SPEC

namespace numpack 
{
namespace SLA
{

//=============================================================================
//Specialization for a Block Sparse Vector
template<class TV>
class SparseVector< DLA::VectorD<TV> > : public SparseLinAlgType< SparseVector< DLA::VectorD<TV> >, false >
{
public:
  typedef DLA::VectorD<TV> Ttype;
  typedef DLA::VectorDView<TV> VectorView_type;

  struct size_type
  {
    size_type(const int m, const int *block_m, const int *block_i) :
      m_(m), block_m_(block_m), block_i_(block_i) {}

    int value_size() const { return block_i_[m_-1] + block_m_[m_-1]; }
    const int m_;
    const int *block_m_;
    const int *block_i_;
  };

  //Used to construct a vector with uniform block sizes
  SparseVector( const int m, const int block_m ) :
    m_(m), block_m_( new int[m] ), block_i_( new int[m] ), values_( new TV[m_*block_m] )
  {
    initialize_block_sizes(m, block_m);
  }

protected:
  void initialize_block_sizes(const int m, const int block_m)
  {
    if (m == 0) return; //zero-length vector

    block_m_[0] = block_m;
    block_i_[0] = 0;
    for (int i = 1; i < m; i++)
    {
      block_m_[i] = block_m;
      block_i_[i] = block_i_[i-1] + block_m_[i];
    }
  }

public:
  // cppcheck-suppress noExplicitConstructor
  SparseVector( const size_type& size ) :
    m_(size.m_), block_m_( new int[m_] ), block_i_( new int[m_] ), values_( new TV[size.value_size()] )
  {
    const int m = m_;
    for (int i = 0; i < m; i++)
    {
      block_m_[i] = size.block_m_[i];
      block_i_[i] = size.block_i_[i];
    }
  }

  SparseVector( const SparseVector& v ) :
    m_(v.m_), block_m_( new int[m_] ), block_i_( new int[m_] ), values_( new TV[v.value_size()] )
  {
    const int m = m_;
    for (int i = 0; i < m; i++)
    {
      block_m_[i] = v.block_m_[i];
      block_i_[i] = v.block_i_[i];
    }
    *this = v;
  }

/* Still need to work this out
  template<class Expr, bool useRF>
  SparseVector( const SparseLinAlgType<Expr, useRF>& r ) :
    m_(v.m_), block_m_( new int[m_] ), block_i_( new int[m_] ), values_( new TV[v.value_size()] )
  {
    *this = r;
  }
*/

  friend class DLA::VectorD< SparseVector< DLA::VectorD<TV> > >;
#ifdef __INTEL_COMPILER
  friend class initializer_list_assign< SparseVector< DLA::VectorD<TV> > >;
#endif
protected:
  SparseVector() : m_(0), block_m_( 0 ), block_i_( 0 ), values_( 0 ) {}
  SparseVector& operator=(const std::initializer_list<int>& s )
  {
    if (s.size() == 2)
      resize(*s.begin(), *(s.begin()+1));
    else
      SANS_ASSERT(false);
    return *this;
  }

  void deallocate()
  {
    delete [] values_;  values_  = NULL;
    delete [] block_m_; block_m_ = NULL;
    delete [] block_i_; block_i_ = NULL;
  }

public:
  ~SparseVector()
  {
    deallocate();
  }

  //Changes the size of the vector by removing all the previous data
  void resize( const int m, const int block_m )
  {
    deallocate();

    m_ = m;
    block_m_ = new int[m];
    block_i_ = new int[m];
    values_ = new TV[m_*block_m];

    initialize_block_sizes(m, block_m);
  }

  void resize( const size_type& size )
  {
    deallocate();

    m_ = size.m_;
    block_m_ = new int[m_];
    block_i_ = new int[m_];
    values_ = new TV[size.value_size()];

    for (int i = 0; i < m_; i++)
    {
      block_m_[i] = size.block_m_[i];
      block_i_[i] = size.block_i_[i];
    }
  }

  int block_m( const int& i ) const { return block_m_[i]; }
  int block_i( const int& i ) const { return block_i_[i]; }
  int m()                     const { return m_; }
  int n()                     const { return 1; }
  size_type size()            const { return size_type(m_, block_m_, block_i_); }
  int value_size()            const { return block_i_[m_-1] + block_m_[m_-1]; }

  SparseVector& operator=(const Real v);
  SparseVector& operator=(const SparseVector& v);
  SparseVector& operator/=(const Real v);
  SparseVector& operator*=(const Real v);

  //L2 norm
  Real norm() const
  {
    Real norm = 0;
    for (int i = 0; i < m_; i++)
      norm += numpack::dot( (*this)[i],(*this)[i] );

    return sqrt(norm);
  }

  //Dot product with another vector
  Real dot( const SparseVector& v ) const
  {
    SANS_ASSERT( m_ == v.m() );

    Real dot = 0;
    for (int i = 0; i < m_; i++)
      dot += numpack::dot( (*this)[i], v[i] );

    return dot;
  }

  //Access operators
        VectorView_type operator[](const unsigned int i)       { return VectorView_type(values_ + block_i_[i], block_m_[i], 1); }
  const VectorView_type operator[](const unsigned int i) const { return VectorView_type(values_ + block_i_[i], block_m_[i], 1); }

  //Lazy evaluation functions
  void value(const Real sgn, SparseVector& b) const;
  void plus(const Real sgn, SparseVector& b) const;

  //Lazy assignment operators
  template< class Expr > inline SparseVector& operator=( const SparseLinAlgType<Expr, true>& r );
  template< class Expr > inline SparseVector& operator+=( const SparseLinAlgType<Expr, true>& r );
  template< class Expr > inline SparseVector& operator-=( const SparseLinAlgType<Expr, true>& r );

  template< class Expr > inline SparseVector& operator=( const SparseLinAlgType<Expr, false>& r );
  template< class Expr > inline SparseVector& operator+=( const SparseLinAlgType<Expr, false>& r );
  template< class Expr > inline SparseVector& operator-=( const SparseLinAlgType<Expr, false>& r );

protected:
  int m_;        //Number of rows in the vector
  int *block_m_; //Number of rows in each of the 'block' vectors
  int *block_i_; //Index in memory for each 'block' vector
  TV *values_;   //Values of the vector
};

} //namespace SLA
} //namespace numpack 

#define TV_SPEC DLA::VectorD<TV>
#include "SparseVector_impl.h"
#undef TV_SPEC

namespace numpack 
{
namespace SLA
{

// I/O
template <class T>
std::ostream&
operator<<( std::ostream& out, const SparseVector<T>& v )
{
  for (int i = 0; i < v.m()-1; i++)
    out << v[i] << ", " << std::endl;
  if (v.m() > 0 )
    out << v[v.m()-1];
  return out;
}

} //namespace SLA
} //namespace numpack 


#endif //SPARSEVECTOR_SCALAR_H
