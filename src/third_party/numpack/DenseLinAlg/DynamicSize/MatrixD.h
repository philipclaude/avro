// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXD_CLASS_H
#define MATRIXD_CLASS_H

#include <iostream>
#include <stdint.h> // uintptr_t

// Surreals are registered with boost::is_arithmetic.
// not sure how to register them with std::is_arithmetic
//#include <boost/type_traits/is_arithmetic.hpp>
#include <type_traits>


#include "tools/SANSnumerics.h"     // Real
#include "tools/SANSException.h"
#include "tools/SANSTraitsInitListAssign.h"
#include "tools/SANSTraitsScalar.h"
#include "tools/init_array_copy.h"

#include "numpack/MatrixScatterAdd.h"

#include "numpack/DenseLinAlg/tools/Identity.h"
#include "numpack/DenseLinAlg/tools/Matrix_Util.h"

#include "MatrixD_Type.h"

#include "MatrixD_Mul.h"
#include "MatrixD_Add.h"
#include "MatrixD_NonZeroPattern.h"

#include "MatrixD_TupleExpression.h"
#include "MatrixD_TupleMatrix.h"

namespace numpack
{

namespace SLA
{
//Forward declare
template<class INT, class T>
class ScalarMatrix_CRS;
}

#ifdef __INTEL_COMPILER
//
// This is completely unnessary if the intel compiler could use templated initializer_list functions....
//
template<class T>
struct initializer_list_assign< DLA::MatrixD< T > >
{
  initializer_list_assign(DLA::MatrixD< T >& val, const std::initializer_list<int>& s ) { val = s; }
};
#endif

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
class MatrixDView : public MatrixDType< MatrixDView<T>, !std::is_arithmetic<T>::value >,
                    public MatrixScatterAdd<T>
{
protected:
  // Needed create MatrixD< MatrixD<T> >
  MatrixDView() : MatrixScatterAdd<T>(LA::eMatrixD, 0, 0), v_ (NULL), stride_(0) {}

public:
  typedef T node_type;

  struct size_type
  {
    size_type(const T* v, const int m, const int n) : v_(v), m_(m), n_(n) {}
    const T* v_;
    const int m_;
    const int n_;
    operator int() const { return m_*n_; }
  };

  explicit MatrixDView( T v0[], const int m, const int n ) :
    MatrixScatterAdd<T>(LA::eMatrixD, m, n), v_(v0), stride_(n)
  {
    SANS_ASSERT( m_ >= 0 );
    SANS_ASSERT( n_ >= 0 );
  }
  explicit MatrixDView( T v0[], const int m, const int n, const int stride ) :
    MatrixScatterAdd<T>(LA::eMatrixD, m, n), v_(v0), stride_(stride)
  {
    SANS_ASSERT( m_ >= 0 );
    SANS_ASSERT( n_ >= 0 );
  }
  MatrixDView( T v0[], const int m, const int n, const typename POD<T>::type& s ) :
    MatrixScatterAdd<T>(LA::eMatrixD, m, n), v_(v0), stride_(n)
  {
    SANS_ASSERT( m_ >= 0 );
    SANS_ASSERT( n_ >= 0 );
    operator=(s);
  }
  MatrixDView( const MatrixDView& z ) : MatrixScatterAdd<T>(z), v_(z.v_), stride_(z.stride_) {}
  ~MatrixDView() {}

  // post creation initialization mainly for arrays of MatrixDView
  //void init( T v0[], const int m, const int n ) { v_ = v0; m_ = m; n_ = n; stride_ = n; SANS_ASSERT( m_ > 0 ); SANS_ASSERT( n_ > 0 ); }

  // matrix dimensions and stride of the array
  size_type size() const { return size_type(v_,m_,n_); }
  using MatrixScatterAdd<T>::m;
  using MatrixScatterAdd<T>::n;
  int stride() const { return stride_; }
        T* data()       { return v_; }
  const T* data() const { return v_; }

  // value accessor operators
        T& operator()( const int i, const int j )       { return v_[i*stride_ + j]; }
  const T& operator()( const int i, const int j ) const { return v_[i*stride_ + j]; }

  //Const versions
  MatrixDView row( const int i ) const { return MatrixDView(v_ + i*stride_,  1, n_, stride_ ); }
  MatrixDView col( const int j ) const { return MatrixDView(v_ + j        , m_,  1, stride_ ); }
  MatrixDView sub( const int i, const int j, const int m, const int n ) const
  {
    return MatrixDView(v_ + i*stride_ + j, m, n, stride_ );
  }
  VectorDView<T> subcol( const int i, const int j, const int n ) const;

  //non-Const versions
  MatrixDView row( const int i ) { return MatrixDView(v_ + i*stride_,  1, n_, stride_ ); }
  MatrixDView col( const int j ) { return MatrixDView(v_ + j        , m_,  1, stride_ ); }
  MatrixDView sub( const int i, const int j, const int m, const int n )
  {
    return MatrixDView(v_ + i*stride_ + j, m, n, stride_ );
  }
  VectorDView<T> subcol( const int i, const int j, const int n );

  // assignment
  MatrixDView& operator=( const MatrixDView& );
  MatrixDView& operator=( const T& );
  MatrixDView& operator=( const typename POD<T>::type& );
  MatrixDView& operator=( const Identity& );
  MatrixDView& operator=( const std::initializer_list< std::initializer_list<T> >& s );

  // Assignment and binary accumulation operators with recursive function expression templates
  template<class Expr> MatrixDView& operator= ( const MatrixDType<Expr, true>& );
  template<class Expr> MatrixDView& operator+=( const MatrixDType<Expr, true>& );
  template<class Expr> MatrixDView& operator-=( const MatrixDType<Expr, true>& );
  template<class Expr, bool useRF> MatrixDView& operator*=( const MatrixDType<Expr, useRF>& );

  // Assignment and binary accumulation operators with element-wise expression templates
  template<class Expr> MatrixDView& operator= ( const MatrixDType<Expr, false>& );
  template<class Expr> MatrixDView& operator+=( const MatrixDType<Expr, false>& );
  template<class Expr> MatrixDView& operator-=( const MatrixDType<Expr, false>& );

  // unary operators; no side effects
  const MatrixDView& operator+() const;

  // binary accumulation operators
  MatrixDView& operator*=( const T& );
  //MatrixDView& operator/=( const T& );

  MatrixDView& operator*=( const typename POD<T>::type& );
  MatrixDView& operator/=( const typename Scalar<T>::type& );

  // input/output
  template<class> friend std::istream& operator>>( std::istream&, MatrixDView<T>& );
  template<class> friend std::ostream& operator<<( std::ostream&, const MatrixDView<T>& );

  //Row manipulations
  inline void swap_rows(const int i, const int j );
  template<class Ta>
  inline void scale_row(const int i, const Ta& a, const int start = 0, int end = 0);
  template<class Ta>
  inline void axpy_rows(const int i, const int j, const Ta& a, const int start = 0, int end = 0);
  inline int max_row_in_col(const int j, const int start = 0) const;

  // Sparse matrix add to the matrix
  void scatterAdd( const MatrixDView& M, const int Map[], int nMap );
  void scatterAdd( const MatrixDView& M, const int rowMap[], int nRow, const int colMap[], int nCol );

  //Allow DenesMatrix to access members directly
  friend class MatrixD<T>;

  //Lazy Expression Operators
  // Assign the value of this matrix to res
  inline void value(const Real sgn, MatrixDView& res) const
  {
    if ( this == &res && sgn == 1 ) return; //Assigning my this to this, nothing to do.

    SANS_ASSERT( this != &res ); //The variable on the left of the assignment cannot also appear on the right
//    SANS_ASSERT( m_ == res.m_ );
//    SANS_ASSERT( n_ == res.n_ );

    for (int i = 0; i < m_; ++i)
      for (int j = 0; j < n_; ++j)
        res(i,j) = sgn*(*this)(i,j);
  }

  inline void plus(const Real sgn, MatrixDView& res) const
  {
    SANS_ASSERT( this != &res ); //The variable on the left of the assignment cannot also appear on the right
//    SANS_ASSERT( m_ == res.m_ );
//    SANS_ASSERT( n_ == res.n_ );

    for (int i = 0; i < m_; ++i)
      for (int j = 0; j < n_; ++j)
        res(i,j) += sgn*(*this)(i,j);
  }

  //Lazy Expression Operators to work with matrix tuples
  // Assign the value of this matrix to res
  template< class MatrixL >
  inline void value(const Real sgn, MatrixDTuple<MatrixL>& res) const
  {
//    SANS_ASSERT( m_ == res.m() );
//    SANS_ASSERT( n_ == res.n() );

    MatrixDView<T> LM( const_cast<T*>(&(*this)(0,       0)), res.Lm_, res.Ln_, stride_ );
    MatrixDView<T> RM( const_cast<T*>(&(*this)(0, res.Ln_)), res.Rm_, res.Rn_, stride_ );
    res.L = sgn*LM;
    res.R = sgn*RM;
  }
  template< class MatrixL >
  inline void plus(const Real sgn, MatrixDTuple<MatrixL>& res) const
  {
//    SANS_ASSERT( m_ == res.m() );
//    SANS_ASSERT( n_ == res.n() );

    MatrixDView<T> LM( const_cast<T*>(&(*this)(0,       0)), res.Lm_, res.Ln_, stride_ );
    MatrixDView<T> RM( const_cast<T*>(&(*this)(0, res.Ln_)), res.Rm_, res.Rn_, stride_ );
    res.L += sgn*LM;
    res.R += sgn*RM;
  }

  //Element wise expression
  inline const node_type& value(const int i) const
  {
    return v_[i];
  }

  //A unique ID for the matrix so the matrix can be uniquely identified
  uintptr_t ID() const { return (uintptr_t)v_; }

  void dump( int indentSize=0, std::ostream& out = std::cout ) const;

protected:
  T *v_;       // value
  using MatrixScatterAdd<T>::m_; // m X n dimensions of the matrix
  using MatrixScatterAdd<T>::n_;
  int stride_; // Stride for when the memory is larger than m x n
};

//----------------------------------------------------------------------------//
// MatrixD: Linear Algebra Dense Matrix
//
// MatrixD allocates/deallocates memory
//
//----------------------------------------------------------------------------//
template< class T >
class MatrixD : public MatrixDView< T >
//class MatrixD : public MatrixDType< MatrixD<T>, true >
{
public:
  typedef T node_type;
  typedef DenseNonZeroPattern<T> NonZeroPattern;

  //Constructors

  friend class MatrixD< MatrixD<T> >;
  friend class MatrixD< DenseNonZeroPattern<T> >;
  typedef typename MatrixDView< T >::size_type size_type;

  // Needed create MatrixD< MatrixD<T> >
  MatrixD() : MatrixDView<T>() {}

protected:

  //Used in DLA::MatrixD to allocate an array of MatrixD
  MatrixD& operator=(const DenseNonZeroPattern<T>& nz)
  {
    delete [] v_; v_ = NULL;
    v_ = new T[nz.nnz()];
    m_ = nz.m();
    n_ = nz.n();
    stride_ = n_;
    return *this;
  }

public:
  //Don't initialize any data here. This treats the matrix more like a POD data type
  //and allows valgrind to catch use of uninitialized data
  explicit MatrixD( const int m, const int n ) : MatrixDView<T>(new T[m*n], m, n) {}
  explicit MatrixD( const DenseMatrixSize& s ) : MatrixDView<T>(new T[s.m()*s.n()], s.m(), s.n()) {}
  explicit MatrixD( const DenseNonZeroPattern<T>& nz ) : MatrixDView<T>(new T[nz.nnz()], nz.m(), nz.n()) {}
  explicit MatrixD( const size_type& s ) { resize(s); }

  MatrixD( const MatrixDView<T>& M ) : MatrixDView<T>(new T[M.m_*M.n_], M.m_, M.n_)
  {
    init_array_copy(m_*n_, M.v_, v_);
    for (int i = 0; i < m_; ++i)
      for (int j = 0; j < n_; ++j)
        (*this)(i,j) = M(i,j);
  }

  MatrixD( const MatrixD& M ) : MatrixDView<T>(new T[M.m_*M.n_], M.m_, M.n_)
  {
    init_array_copy(m_*n_, M.v_, v_);
    for (int i = 0; i < m_; ++i)
      for (int j = 0; j < n_; ++j)
        (*this)(i,j) = M(i,j);
  }

  // move constructor just takes the memory
  MatrixD( MatrixD&& M ) : MatrixDView<T>(M.v_, M.m_, M.n_)
  {
    M.v_ = NULL;
  }

  template<class T2>
  MatrixD( const MatrixDView<T2>& M ) : MatrixDView<T>(new T[M.m()*M.n()], M.m(), M.n())
  {
    init_array_copy(m_*n_, &M(0,0), v_);
    for (int i = 0; i < m_; ++i)
      for (int j = 0; j < n_; ++j)
        (*this)(i,j) = M(i,j);
  }

  explicit MatrixD( const int m, const int n, const T v0[] ) : MatrixDView<T>(new T[m*n], m, n)
  {
    init_array_copy(m_*n_, v0, v_);
    const int size = m*n;
    for (int i = 0; i < size; i++)
      v_[i] = v0[i];
  }

  explicit MatrixD( const int m, const int n, const T& v0 ) : MatrixDView<T>(new T[m*n], m, n)
  {
    const int size = m*n;
    for (int i = 0; i < size; i++)
      v_[i] = v0;
  }

  explicit MatrixD( const int m, const int n, const typename POD<T>::type& v0 ) : MatrixDView<T>(new T[m*n], m, n)
  {
    const int size = m*n;
    for (int i = 0; i < size; i++)
      v_[i] = v0;
  }

  MatrixD( const std::initializer_list< std::initializer_list<T> >& s )
    : MatrixDView<T>(new T[s.size()*(*s.begin()).size()], s.size(), (*s.begin()).size())
  {
    int k = 0;
    for (auto i = s.begin(); i != s.end(); ++i)
    {
      SANS_ASSERT( (int)(*i).size() == n_ );
      init_array_copy(n_, (*i).begin(), v_+k);
      for (auto j = (*i).begin(); j != (*i).end(); ++j)
        v_[k++] = *j;
    }
  }

  template<class U>
  MatrixD( const std::initializer_list< std::initializer_list<U> >& s )
    : MatrixDView<T>(new T[s.size()*(*s.begin()).size()], s.size(), (*s.begin()).size())
  {
    int k = 0;
    for (auto i = s.begin(); i != s.end(); ++i)
    {
      SANS_ASSERT( (int)(*i).size() == n_ );
      init_array_copy(n_, (*i).begin(), v_+k);
      for (auto j = (*i).begin(); j != (*i).end(); ++j)
        v_[k++] = *j;
    }
  }

#if __INTEL_COMPILER
  //Maybe someday the intel compiler will get fixed and we won't need any of this mess...
  MatrixD( const std::initializer_list< std::initializer_list< std::initializer_list<int> > >& s )
    : MatrixDView<T>(new T[s.size()*(*s.begin()).size()], s.size(), (*s.begin()).size())
  {
    int k = 0;
    for (auto i = s.begin(); i != s.end(); ++i)
    {
      SANS_ASSERT( (int)(*i).size() == n_ );
      for (auto j = (*i).begin(); j != (*i).end(); ++j)
      {
        initializer_list_assign<T>(v_[k], *j);
        k++;
      }
    }
  }
#else
  template<class U>
  MatrixD( const std::initializer_list< std::initializer_list< std::initializer_list<U> > >& s )
    : MatrixDView<T>(new T[s.size()*(*s.begin()).size()], s.size(), (*s.begin()).size())
  {
    int k = 0;
    for (auto i = s.begin(); i != s.end(); ++i)
    {
      SANS_ASSERT( (int)(*i).size() == n_ );
      init_array_copy(n_, (*i).begin(), v_+k);
      for (auto j = (*i).begin(); j != (*i).end(); ++j)
        v_[k++] = *j;
    }
  }
#endif

  //Constructs a dense matrix from a sparse matrix
  MatrixD( const SLA::ScalarMatrix_CRS<int, Real>& A );

  //This allows a matrix to be assigned an expression as it is constructed, i.e.
  //MatrixD<T> C = A + B;
  template<class Expr, bool useRF>
  MatrixD( const MatrixDType<Expr, useRF>& r )
    : MatrixDView<T>(new T[r.cast().size()],
      r.cast().m(),
      r.cast().n())
  {
    this->operator=(r);
  }

  ~MatrixD() { delete [] v_; }

  // assignment operator. All other methods are inherited
  MatrixD& operator=( const MatrixD& M )
  {
    // Allow the matrix to be resized if it's NULL
    if (v_ == NULL) resize(M.size());

    SANS_ASSERT( m_ == M.m_ );
    SANS_ASSERT( n_ == M.n_ );

    for (int i = 0; i < m_; ++i)
      for (int j = 0; j < n_; ++j)
        (*this)(i,j) = M(i,j);

    return *this;
  }

  template<class Expr, bool useRF>
  MatrixD& operator=( const MatrixDType<Expr, useRF>& r )
  {
    // Allow the matrix to be resized if it's NULL
    if (v_ == NULL) resize(r.cast().m(), r.cast().n());

    MatrixDView<T>::operator=(r);
    return *this;
  }

  MatrixD& operator=( MatrixD&& M )
  {
    if (v_ == NULL)
    {
      // Simply steal the memory
      m_ = M.m();
      n_ = M.n();
      stride_ = n_;
    }
    else
    {
      // Check the size match before stealing the memory
      SANS_ASSERT( m_ == M.m_ );
      SANS_ASSERT( n_ == M.n_ );
      delete [] v_; v_ = NULL;
    }

    // Steal the memory
    v_ = M.v_;
    M.v_ = NULL;

    return *this;
  }


  MatrixD& operator=( const T& v0 ) { MatrixDView<T>::operator=(v0); return *this; }
  MatrixD& operator=( const typename POD<T>::type& v0 ) { MatrixDView<T>::operator=(v0); return *this; }
  MatrixD& operator=( const Identity& I ) { MatrixDView<T>::operator=(I); return *this; }

  //Special assignment operator used to resize MatrixD<MatrixD>
  MatrixD& operator=( const std::initializer_list<int>& s )
  {
    SANS_ASSERT(v_ == NULL);
    SANS_ASSERT(m_ == 0);
    SANS_ASSERT(n_ == 0);
    SANS_ASSERT(s.size() == 2);
    m_ = *s.begin();
    n_ = *(s.begin()+1);
    SANS_ASSERT(m_ >= 0);
    SANS_ASSERT(n_ >= 0);
    stride_ = n_;

    v_ = new T[m_*n_];

    return *this;
  }

  void resize( const size_type& s )
  {
    delete [] v_; v_ = NULL;
    v_ = new T[s.m_*s.n_];
    m_ = s.m_;
    n_ = s.n_;
    stride_ = n_;
    init_array_copy(s.m_*s.n_, s.v_, v_);
  }

  void resize( const int m, const int n )
  {
    delete [] v_; v_ = NULL;
    v_ = new T[m*n];
    m_ = m;
    n_ = n;
    stride_ = n_;
  }

protected:
  using MatrixDView<T>::v_;
  using MatrixDView<T>::m_;
  using MatrixDView<T>::n_;
  using MatrixDView<T>::stride_;
};


// assignment

template< class T >
inline MatrixDView<T>&
MatrixDView<T>::operator=( const MatrixDView<T>& M )
{
  SANS_ASSERT_MSG( m_ == M.m_, "%d == %d", m_, M.m_ );
  SANS_ASSERT_MSG( n_ == M.n_, "%d == %d", n_, M.n_ );

  for (int i = 0; i < m_; i++)
    for (int j = 0; j < n_; j++)
      (*this)(i,j) = M(i,j);

  return *this;
}

template< class T >
inline MatrixDView<T>&
MatrixDView<T>::operator=( const T& r )
{
  for (int i = 0; i < m_; i++)
    for (int j = 0; j < n_; j++)
      (*this)(i,j) = r;

  return *this;
}

//This is needed when T is not POD. It will define an operator= for POD
template< class T >
inline MatrixDView<T>&
MatrixDView<T>::operator=( const typename POD<T>::type& r )
{
  for (int i = 0; i < m_; i++)
    for (int j = 0; j < n_; j++)
      (*this)(i,j) = r;

  return *this;
}

template< class T >
inline MatrixDView<T>&
MatrixDView<T>::operator=( const Identity& I )
{
  *this = 0;
  if ( m_ == n_ )
    for (int i = 0; i < m_; i++)
      (*this)(i,i) = I;
  else
    for (int i = 0; i < m_; i++)
      for (int j = 0; j < n_; j++)
        if ( i == j )
          (*this)(i,j) = I;

  return *this;
}

template< class T >
inline MatrixDView<T>&
MatrixDView<T>::operator=( const std::initializer_list< std::initializer_list<T> >& s )
{
  SANS_ASSERT( (int)s.size() == m_ );
  int ii = 0;
  for (typename std::initializer_list< std::initializer_list<T> >::iterator i = s.begin(); i != s.end(); ++i)
  {
    SANS_ASSERT( (int)(*i).size() == n_ );
    int jj = 0;
    for (typename std::initializer_list<T>::iterator j = (*i).begin(); j != (*i).end(); ++j)
      (*this)(ii,jj++) = *j;
    ii++;
  }

  return *this;
}


// Assignment and binary accumulation operators with recursive function expression templates

template< class T >
template< class Expr >
inline MatrixDView<T>&
MatrixDView<T>::operator=( const MatrixDType<Expr, true>& r )
{
  const Expr& Tree = r.cast();

  SANS_ASSERT( m_ == Tree.m() );
  SANS_ASSERT( n_ == Tree.n() );

  Tree.value(1, *this);

  return *this;
}

template< class T >
template< class Expr >
inline MatrixDView<T>&
MatrixDView<T>::operator+=( const MatrixDType<Expr, true>& r )
{
  const Expr& Tree = r.cast();

  SANS_ASSERT( m_ == Tree.m() );
  SANS_ASSERT( n_ == Tree.n() );

  Tree.plus(1, *this);

  return *this;
}

template< class T >
template< class Expr >
inline MatrixDView<T>&
MatrixDView<T>::operator-=( const MatrixDType<Expr, true>& r )
{
  const Expr& Tree = r.cast();

  SANS_ASSERT( m_ == Tree.m() );
  SANS_ASSERT( n_ == Tree.n() );

  Tree.plus(-1, *this);

  return *this;
}

template< class T >
template< class Expr, bool useRF >
inline MatrixDView<T>&
MatrixDView<T>::operator*=( const MatrixDType<Expr, useRF>& r )
{
  MatrixD<T> tmp(*this);
  *this = tmp*r;
  return *this;
}


// Assignment and binary accumulation operators with element-wise expression templates

template< class T >
template< class Expr >
inline MatrixDView<T>&
MatrixDView<T>::operator=( const MatrixDType<Expr, false>& r )
{
  const Expr& Tree = r.cast();

  SANS_ASSERT_MSG( m_ == Tree.m() && n_ == Tree.n(),
                  "m_ = %d, n_ = %d, Tree.m() = %d, Tree.n() = %d", m_, n_, Tree.m(), Tree.n() );

  //Saving to a local temporary improves optimization as the compiler does not
  //need to assume that m_ or n_ could change in the loop
  const int tmp = m_*n_;

  for ( int i = 0; i < tmp; i++)
    v_[i] = Tree.value(i);

  return *this;
}

template< class T >
template< class Expr >
inline MatrixDView<T>&
MatrixDView<T>::operator+=( const MatrixDType<Expr, false>& r )
{
  const Expr& Tree = r.cast();

  SANS_ASSERT_MSG( m_ == Tree.m() && n_ == Tree.n(),
                  "m_ = %d, n_ = %d, Tree.m() = %d, Tree.n() = %d", m_, n_, Tree.m(), Tree.n() );

  //Saving to a local temporary improves optimization as the compiler does not
  //need to assume that m_ or n_ could change in the loop
  const int tmp = m_*n_;

  for ( int i = 0; i < tmp; i++)
    v_[i] += Tree.value(i);

  return *this;
}

template< class T >
template< class Expr >
inline MatrixDView<T>&
MatrixDView<T>::operator-=( const MatrixDType<Expr, false>& r )
{
  const Expr& Tree = r.cast();

  SANS_ASSERT_MSG( m_ == Tree.m() && n_ == Tree.n(),
                  "m_ = %d, n_ = %d, Tree.m() = %d, Tree.n() = %d", m_, n_, Tree.m(), Tree.n() );

  //Saving to a local temporary improves optimization as the compiler does not
  //need to assume that m_ or n_ could change in the loop
  const int tmp = m_*n_;

  for ( int i = 0; i < tmp; i++)
    v_[i] -= Tree.value(i);

  return *this;
}

// unary operators; no side effects

template< class T >
inline const MatrixDView<T>&
MatrixDView<T>::operator+() const
{
  return *this;
}

template< class Expr, bool useRF >
inline const OpMulDScalar<Expr, useRF>
operator-(MatrixDType<Expr, useRF> const& e)
{
  return OpMulDScalar<Expr, useRF>( e.cast(), -1 );
}

// binary accumulation operators

template< class T >
inline MatrixDView<T>&
MatrixDView<T>::operator*=( const T& r )
{
  for (int i = 0; i < m_; i++)
    for (int j = 0; j < n_; j++)
      (*this)(i,j) *= r;

  return *this;
}

//This is needed when T is not POD. It will define an operator*= for POD
template< class T >
inline MatrixDView<T>&
MatrixDView<T>::operator*=( const typename POD<T>::type& r )
{
  for (int i = 0; i < m_; i++)
    for (int j = 0; j < n_; j++)
      (*this)(i,j) *= r;

  return *this;
}

#if 0
//Forward declaration
struct InverseLU;

template< class T >
inline MatrixDView<T>&
MatrixDView<T>::operator/=( const T& r )
{
  T tmp = InverseLU::Inverse( r );
  for (int i = 0; i < m_; i++)
    for (int j = 0; j < n_; j++)
      (*this)(i,j) *= tmp;

  return *this;
}
#endif

//This is needed when T is not POD. It will define an operator/= for POD
template< class T >
inline MatrixDView<T>&
MatrixDView<T>::operator/=( const typename Scalar<T>::type& r )
{
  for (int i = 0; i < m_; i++)
    for (int j = 0; j < n_; j++)
      (*this)(i,j) /= r;

  return *this;
}

//-----------------------------------------------------------------------------
    //Row operations

//Swaps rows i and j
template< class T >
inline void
MatrixDView<T>::swap_rows(const int i, const int j )
{
  if ( i == j ) return; //Nothing to do...
  MatrixUtil<T,T>::swap(&(*this)(i,0), &(*this)(j,0), n_);
}

//Scales the row i
template< class T >
template<class Ta>
inline void
MatrixDView<T>::scale_row(const int i, const Ta& a, const int start, int end)
{
  if (end == 0) end = n_;
  MatrixUtil<T,Ta>::scal(&(*this)(i,0) + start, a, end - start);
}

//Performs y = a*x + y where x = row(i) and y = row(j)
template< class T >
template< class Ta >
inline void
MatrixDView<T>::axpy_rows(const int i, const int j, const Ta& a, const int start, int end)
{
  if (end == 0) end = n_;
  MatrixUtil<T,Ta>::axpy(&(*this)(i,0) + start, &(*this)(j,0) + start, a, end - start);
}

//Finds the row with maximum value in a column
template< class T >
inline int
MatrixDView<T>::max_row_in_col(const int j, const int start) const
{
  return MatrixUtil<T,T>::max_row_in_col(&(*this)(0,j), m_, stride_, start);
}

template <class T>
inline void
MatrixDView<T>::scatterAdd( const MatrixDView& M, const int Map[], int nMap )
{
  for (int i = 0; i < nMap; i++)
  {
    int iGlobal = Map[i];
    for (int j = 0; j < nMap; j++)
    {
      int jGlobal = Map[j];
      (*this)(iGlobal,jGlobal) += M(i,j);
    }
  }
}

template <class T>
inline void
MatrixDView<T>::scatterAdd( const MatrixDView& M, const int rowMap[], int nRow, const int colMap[], int nCol )
{
  for (int i = 0; i < nRow; i++)
  {
    int iGlobal = rowMap[i];
    for (int j = 0; j < nCol; j++)
    {
      int jGlobal = colMap[j];
      (*this)(iGlobal,jGlobal) += M(i,j);
    }
  }
}


// I/O

// debug dump of private data
template <class T>
void
MatrixDView<T>::dump( int indentSize, std::ostream& out ) const
{
  std::string indent(indentSize, ' ');
  out << indent << "MatrixD<" << m_ << "," << n_ << ",T>:" << std::endl;
#if 1
  out << indent << "  data = ";
  for (int n = 0; n < m_*n_; n++)
    out << v_[n] << " ";
  out << std::endl;
#else     // only works for class T with member function dump()
  for (int n = 0; n < m_*n_; n++)
  {
    out << indent << "  data[" << n << "] = ";
    data[n].dump(indentSize);
  }
#endif
}

template <>
void
MatrixDView<Real>::dump( int indentSize, std::ostream& out ) const;

template <>
void
MatrixDView< MatrixS<1,1,Real> >::dump( int indentSize, std::ostream& out ) const;

template<class T>
inline void dump( const MatrixDView<T>& mtx, int indentSize )
{
  std::string indent(indentSize, ' ');
  int m = mtx.m();
  int n = mtx.n();
  std::cout << indent << "{";
  for (int i = 0; i < m; i++)
  {
    std::cout << " {";
    for (int j = 0; j < n; j++)
    {
      std::cout << " " << mtx(i,j);
    }
    std::cout << " }";
    if (i < m-1)  std::cout << ",";
  }
  std::cout << std::endl;
}

template <class T>
std::ostream&
operator<<( std::ostream& out, const MatrixDView<T>& M )
{
  for (int i = 0; i < M.m(); i++)
  {
    for (int j = 0; j < M.n(); j++)
      out << M(i,j) << ", ";
    if ( i < M.m()-1 ) out << std::endl;
  }
  return out;
}

template <class T>
std::istream&
operator>>( std::istream &in, MatrixDView<T>& M )
{
  for (int i = 0; i < M.m(); i++)
    for (int j = 0; j < M.n(); j++)
      in >> M(i,j);
  return in;
}

} //namespace DLA
} //namespace numpack


#endif // MATRIXD_CLASS_H
