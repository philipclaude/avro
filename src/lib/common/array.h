#ifndef avro_LIB_ARRAY_H_
#define avro_LIB_ARRAY_H_

#include "common/error.h"
#include "common/types.h"

#include <vector>

namespace avro
{

template<typename type>
class Array
{
public:
  Array()
  {}

  type& operator[]       (index_t k) { return data_[k]; }
  const type operator[] (index_t k) const { return data_[k]; }

  index_t nb() const { return data_.size(); }

  void set( index_t k , const type& x )
  {
    avro_assert( k < nb() );
    data_[k] = x;
  }

  void add( type* x , index_t n )
  {
    for (index_t j=0;j<n;j++)
      add(x[j]);
  }

  void add( const type* x , index_t n )
  {
    for (index_t j=0;j<n;j++)
      add(x[j]);
  }

  void add( const type x ) { data_.push_back(x); }

  void insert( index_t loc , const type* x , index_t n )
  {
    data_.insert( data_.begin()+loc , x , x+n );
  }

  void insert( index_t loc , const type x )
  {
    data_.insert( data_.begin()+loc , &x , &x+1 );
  }

  void remove( index_t k0 )
  {
    data_.erase( data_.begin() + k0 );
  }

  std::vector<type> data() const { return data_; }
  std::vector<type>& get_data() const { return data_; }

  void clear()
  {
    data_.clear();
  }

  void resize( index_t n , type x )
  {
    data_.resize( n , x );
  }

  void print() const
  {
    avro_implement;
  }

private:

  void remove( index_t start , index_t end )
  {
    data_.erase( data_.begin() + start , data_.begin() + end );
  }

protected:

  std::vector<type> data_;
};

} // avro

#endif
