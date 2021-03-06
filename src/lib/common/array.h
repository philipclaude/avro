//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_ARRAY_H_
#define avro_LIB_ARRAY_H_

#include "common/error.h"
#include "avro_types.h"

#include <set>
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
    avro_assert_msg( k < nb() , "k = %lu but nb = %lu", k , nb() );
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

  void reserve( index_t n ) {
    data_.reserve(n);
  }

  void shrink_to_fit() {
    data_.shrink_to_fit();
  }

  void batch_erase( index_t n ) {
    std::vector<type> data0(data_.begin()+n,data_.end());
    data_.assign(data0.begin(),data0.end());
  }

  void move_to_front( const std::vector<index_t>& idx ) {

    // make enough space for the data
    std::vector<type> data;
    data.reserve(data_.size());

    // move the requested indices to the front
    for (index_t k = 0; k < idx.size(); k++) {
      data.push_back(data_[idx[k]]);
    }

    // now add the other data behind the stuff we moved to the front
    std::set<index_t> sidx(idx.begin(),idx.end());
    for (index_t k = 0; k < nb(); k++) {
      if (sidx.find(k) != sidx.end()) continue; // this element was moved to the front
      data.push_back(data_[k]);
    }
    avro_assert( data_.size() == data.size() );
    data_.assign( data.begin() , data.end() );
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
