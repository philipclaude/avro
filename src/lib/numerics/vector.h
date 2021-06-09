#ifndef PYDG_NUMERICS_VECTOR_H_
#define PYDG_NUMERICS_VECTOR_H_

#include "common/error.h"
#include "common/types.h"

#include "numerics/matrix.h"

#include <vector>

namespace avro
{

template<typename type>
class vecd
{
public:
  vecd(index_t m) :
    m_(m),
    data_(m)
  {}

  type& operator() (index_t i) {
    avro_assert_msg( i < m_ , "attempt to access i = %lu but m = %lu" , i , m_ );
    return data_[i];
  }

  const type& operator() (index_t i) const {
    avro_assert_msg( i < m_ , "attempt to access i = %lu but m = %lu" , i , m_ );
    return data_[i];
  }

  // zero out the matrix
  void zero() {
    for (index_t i = 0; i < m_; i++)
      data_[i] = 0;
  }

  index_t m() const { return m_; }

  void print() const {
    printf("%s:\n",__PRETTY_FUNCTION__);
    for (index_t i = 0; i < m_; i++)
      std::cout << "(" + std::to_string(i) + "): " << (*this)(i) << std::endl;
  }

private:
  index_t m_;
  std::vector<type> data_;
};

template<index_t _M,typename type>
class vecs : public mats<_M,1,type>
{
public:
  static const index_t M = _M;

private:
  using mats<M,1,type>::data_;

public:
  vecs() {}
  vecs( const type* data , index_t _m ) {
    avro_assert( _m == M );
    for (index_t i = 0; i < M; i++)
      data_[i] = data[i];
  }

  vecs( int a ) {
    for (index_t i = 0; i < M; i++)
      data_[i] = a;
  }

  // assignment
  template<typename S>
  vecs<M,type>& operator= (const vecs<M,S>& b) {
    for (index_t i = 0; i < M; i++)
      (*this)(i) = b(i);
    return *this;
  }

  vecs<M,type>& operator= (int a) {
    for (index_t i = 0; i < M; i++)
      data_[i] = a;
    return *this;
  }

  type& operator() (index_t i) {
    avro_assert( i < M );
    return data_[i];
  }

  const type& operator() (index_t i) const {
    avro_assert( i < M );
    return data_[i];
  }

  void print() const {
    printf("%s:\n",__PRETTY_FUNCTION__);
    for (index_t i = 0; i < M; i++)
      std::cout << "(" + std::to_string(i) + "): " << (*this)(i) << std::endl;
  }
};

template<typename R,typename S,typename T,index_t M>
T dot( const vecs<M,R>& u , const vecs<M,S>& v );

} // pydg

#endif
