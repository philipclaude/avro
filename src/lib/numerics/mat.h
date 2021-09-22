//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef PYDG_NUMERICS_MATRIX_H_
#define PYDG_NUMERICS_MATRIX_H_

#include "common/error.h"
#include "avro_types.h"

#include "numerics/surreal/SurrealS.h"
#include "numerics/result_of.h"

#include <vector>

namespace avro
{

template<typename T> class vecd;
template<typename T> class symd;

template<typename T>
class matd {
public:
  matd(index_t n) :
    m_(n),
    n_(n),
    data_(n*n,0)
  {}

  matd() :
    m_(0),
    n_(0)
  {}

  matd(index_t m,index_t n) :
    m_(m),
    n_(n),
    data_(m*n,0)
  {}

  matd(const matd& A) :
    m_(A.m()),
    n_(A.n()),
    data_(m_*n_) {
    for (index_t i = 0; i < m_; i++)
    for (index_t j = 0; j < n_; j++)
      (*this)(i,j) = A(i,j);
  }

  matd( const T* v , index_t m , index_t n ) :
    m_(m),
    n_(n),
    data_(v,v+m*n)
  {}

  matd( const symd<T>& A ) :
    m_(A.m()),
    n_(A.m()),
    data_(m_*n_) {
    for (index_t i = 0; i < m_; i++)
    for (index_t j = 0; j < n_; j++)
      (*this)(i,j) = A(i,j);
  }

  index_t m() const { return m_; }
  index_t n() const { return n_; }

  // zero out the matrix
  void zero() {
    for (index_t i = 0; i < m_*n_; i++)
      data_[i] = 0;
  }

  // identity
  void eye() {
    avro_assert( m_ == n_ );
    for (index_t i = 0; i < m_; i++)
      (*this)(i,i) = 1;
  }

  void set( const real_t& x ) {
    for (index_t j = 0; j < data_.size(); j++)
      data_[j] = x;
  }

  void set( const matd<T>& A ) {
    avro_assert( A.m() == m_ && A.n() == n_ );
    for (index_t i = 0; i < m_; i++)
    for (index_t j = 0; j < n_; j++)
      (*this)(i,j) = A(i,j);
  }

  void resize( index_t m , index_t n ) {
    m_ = m;
    n_ = n;
    data_.resize(m*n);
  }

  inline T& operator() (index_t i, index_t j) {
    avro_assert_msg( i < m_ && j < n_ , "i = %lu, j = %lu of %lu x %lu matrix" , i,j,m_,n_ );
    return data_[j*m_+i];
  }

  inline const T& operator() (index_t i, index_t j) const {
    avro_assert_msg( i < m_ && j < n_ , "i = %lu, j = %lu of %lu x %lu matrix" , i,j,m_,n_ );
    return data_[j*m_+i];
  }

  matd<T>& operator=( const symd<T>& S ) {
    resize( S.n() , S.n() );
    for (index_t i = 0; i < m_; i++)
    for (index_t j = 0; j < m_; j++)
      (*this)(i,j) = S(i,j);
    return *this;
  }

  void set( const symd<T>& S );

  void set_row( index_t i , const vecd<T>& row );
  void get_row( index_t i , vecd<T>& row ) const;

  void print() const {
    printf("%s:\n",__PRETTY_FUNCTION__);
    for (index_t i = 0; i < m_; i++)
      for (index_t j = 0; j < n_; j++)
        std::cout << "(" + std::to_string(i) + "," + std::to_string(j) + "): " << (*this)(i,j) << std::endl;
  }

  void dump() const { print(); }

  index_t memory() const {
    return m_*n_*sizeof(T);
  }

private:
  index_t m_;
  index_t n_;
  std::vector<T> data_;
};

template<index_t M,index_t N,typename T>
class mats {
public:
  mats() {
    zero();
  }

  // zero out the matrix
  void zero() {
    for (index_t i = 0; i < M*N; i++)
      data_[i] = 0;
  }

  // identity
  void eye() {
    avro_assert( M == N );
    for (index_t i = 0; i < M; i++)
      (*this)(i,i) = 1;
  }

  // assignment
  template<typename S>
  mats<M,N,T>& operator= (const mats<M,N,S>& b) {
    for (index_t i = 0; i < M; i++)
    for (index_t j = 0; j < N; j++)
      (*this)(i,j) = b(i,j);
    return *this;
  }

  mats<M,N,T>& operator= (const real_t& b) {
    for (index_t i = 0; i < M; i++)
    for (index_t j = 0; j < N; j++)
      (*this)(i,j) = b;
    return *this;
  }

  // indexing (non-const)
  T& operator() (index_t i, index_t j) {
    avro_assert( i < M && j < N );
    return data_[j*M+i];
  }

  // indexing (const)
  const T& operator() (index_t i, index_t j) const {
    avro_assert( i < M && j < N );
    return data_[j*M+i];
  }

  void print() const {
    printf("%s:\n",__PRETTY_FUNCTION__);
    for (index_t i = 0; i < M; i++)
      for (index_t j = 0; j < N; j++)
        std::cout << "(" + std::to_string(i) + "," + std::to_string(j) + "): " << (*this)(i,j) << std::endl;
  }

  const T* data() const { return data_; }

protected:
  T data_[M*N];
};

template<index_t M , index_t N, class T>
std::ostream&
operator<<( std::ostream& os, const mats<M,N,T>& z )
{
  z.print();
  return os;
}

template<index_t M, index_t N,class R>
mats<M,N,R>
operator+ ( const mats<M,N,R>& A , const mats<M,N,R>& B ) {
  mats<M,N,R> C; // initializes to zero
  for (index_t i = 0; i < M; i++)
  for (index_t j = 0; j < N; j++)
    C(i,j) = A(i,j) + B(i,j);
  return C;
}

template<typename T> vecd<T> operator* ( const matd<T>& A , const vecd<T>& x );
template<typename R,typename S> matd< typename result_of<R,S>::type > operator* ( const matd<R>& A , const matd<S>& B );
template<typename R,typename S> matd< typename result_of<R,S>::type > operator* ( const matd<R>& A , const S& b );
template<typename R,typename S> matd< typename result_of<R,S>::type > operator* ( const S& a ,  const matd<R>& B );
template<typename R,typename S> matd< typename result_of<R,S>::type > operator+ ( const matd<R>& A , const matd<S>& B );
template<typename R,typename S> matd< typename result_of<R,S>::type > operator- ( const matd<R>& A , const matd<S>& B );

} // avro

#endif
