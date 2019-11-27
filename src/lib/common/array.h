#ifndef LUNA_LIB_ARRAY_H_
#define LUNA_LIB_ARRAY_H_

#include "common/error.h"
#include "common/types.h"

#include <algorithm>
#include <vector>

namespace luna
{

#if 1

enum ArrayLayoutCategory
{
  ArrayLayout_Jagged,
  ArrayLayout_Rectangular,
  ArrayLayout_Simple,
  ArrayLayout_None
};

template<typename type>
class Array
{
private:
  typedef index_t (Array<type>::*layout_ptr)(index_t k) const;

  index_t layout_jagged( index_t k ) const { return first_[k]; }
  index_t layout_rectangular( index_t k ) const { return rank_*k; }
  index_t layout_simple( index_t k ) const { return k; }

  const layout_ptr
  get_layout( ArrayLayoutCategory category ) const
  {
    if (category==ArrayLayout_Jagged)
      return &Array<type>::layout_jagged;
    if (category==ArrayLayout_Rectangular)
      return &Array<type>::layout_rectangular;
    if (category==ArrayLayout_Simple)
      return &Array<type>::layout_simple;
    return NULL;
  }

public:
  Array( ArrayLayoutCategory category ) :
    category_(category),
    sorted_(false),
    rank_(0),
    layout_ptr_(get_layout(category_))
  {}

  Array( const ArrayLayoutCategory& category , index_t rank ) :
    category_(category),
    sorted_(false),
    rank_(rank),
    layout_ptr_(get_layout(category_))
  {}

  void set_sorted( bool sorted ) { sorted_ = sorted; }

  void set_rank( index_t rank )
  {
    luna_assert( category_ == ArrayLayout_Rectangular );
    rank_ = rank;
  }

  index_t rank() const
  {
    luna_assert( category_==ArrayLayout_Rectangular || category_==ArrayLayout_Simple );
    return rank_;
  }

  type* operator[]       (index_t k)
    { return &data_[location(k)]; }
  const type* operator[] (index_t k) const
    { return &data_[location(k)]; }

  type* operator()       (index_t k)
    { return &data_[location(k)]; }
  const type* operator() (index_t k) const
    { return &data_[location(k)]; }

  type& operator()       (index_t k, index_t j)
    { return data_[location(k)+j]; }
  const type operator() (index_t k, index_t j) const
    { return data_[location(k)+j]; }

  index_t nv( index_t k ) const
  {
    if (category_==ArrayLayout_Rectangular) return rank_;
    else if (category_==ArrayLayout_Simple) return 1;
    else if (category_==ArrayLayout_Jagged) return last_[k]-first_[k];
    else luna_assert_not_reached;
  }

  index_t nb() const
  {
    if (category_==ArrayLayout_Rectangular) return data_.size()/rank_;
    else if (category_==ArrayLayout_Simple) return data_.size();
    else if (category_==ArrayLayout_Jagged) return first_.size();
    else luna_assert_not_reached;
    return 0;
  }

  std::vector<type>
  get( index_t k ) const
  {
    std::vector<type> x( nv(k) );
    for (index_t j=0;j<x.size();j++)
      x[j] = data_[location(k)+j];
    return x;
  }

  void set( index_t k , const type* x )
  {
    for (coord_t j=0;j<nv(k);j++)
      data_[location(k)+j] = x[j];
  }

  void set( index_t k , const type& x )
  {
    luna_assert( category_==ArrayLayout_Simple );
    data_[k] = x;
  }

  void add( type* x , index_t n )
  {
    if (category_==ArrayLayout_Rectangular)
      luna_assert( n == rank_ );
    if (category_==ArrayLayout_Jagged)
      first_.push_back( data_.size() );
    for (index_t j=0;j<n;j++)
      data_.push_back(x[j]);
    if (category_==ArrayLayout_Jagged)
      last_.push_back( data_.size() );
  }

  void add( const type* x , index_t n )
  {
    if (category_==ArrayLayout_Jagged)
      first_.push_back( data_.size() );
    for (index_t j=0;j<n;j++)
      data_.push_back(x[j]);
    if (category_==ArrayLayout_Jagged)
      last_.push_back( data_.size() );
  }

  void add( const type x )
  {
    luna_assert( category_ == ArrayLayout_Simple );
    data_.push_back(x);
  }

  void remove( index_t k )
  {
    luna_implement;
  }

  void insert( index_t loc , const type* x , index_t n )
  {
    data_.insert( data_.begin()+loc , x , x+n );
  }

  void insert( index_t loc , const type x )
  {
    data_.insert( data_.begin()+loc , &x , &x+1 );
  }

  std::vector<type> data() const { return data_; }

  void clear()
  {
    first_.clear();
    last_.clear();
    data_.clear();
    sorted_   = false;
    rank_     = 0;
  }

  void print() const
  {
    std::vector<type> x;
    for (index_t k=0;k<nb();k++)
    {
      x = get(k);
      print_inline(x);
    }
  }

private:


  index_t location( index_t k ) const
  {
    if (category_==ArrayLayout_Rectangular) return k*rank_;
    else if (category_==ArrayLayout_Simple) return k;
    else if (category_==ArrayLayout_Jagged) return first_[k];
    else luna_assert_not_reached;
  }

  const ArrayLayoutCategory category_;
  bool sorted_;
  index_t rank_;
  const layout_ptr layout_ptr_;

  std::vector<type>    data_;
  std::vector<index_t> first_;
  std::vector<index_t> last_;

};

#else

template<typename type>
class ArrayBase
{
public:
  void set_sorted( bool x ) { sorted_ = x; }

  void add( const type* x , index_t n )
  {
    assert( !sorted_ );
    for (index_t i=0;i<n;i++)
      data_.push_back( x[i] );
  }

  void add( type* x , index_t n )
  {
    if (sorted_) std::sort( x , x+n );
    for (index_t i=0;i<n;i++)
      data_.push_back( x[i] );
  }

  void insert( index_t loc , const type* x , index_t n )
  {
    data_.insert( data_.begin()+loc , x , x+n );
  }

  const std::vector<type>& data() const { return data_; }

protected:
  ArrayBase( bool sorted=false ) :
    sorted_(sorted)
  {}

  void remove( index_t start , index_t end )
  {
    data_.erase( data_.begin() + start , data_.begin() + end );
  }

  void clear()
  {
    data_.clear();
    sorted_ = false;
  }

  std::vector<type> data_;

private:
  bool sorted_;
};

class Jagged;
class Rectangular;
class Simple;

template<typename Layout,typename type> class Array;

template<typename type>
class Array<Jagged,type> : public ArrayBase<type>
{
public:
  void add( const type* x , index_t n )
  {
    first_.push_back( ArrayBase<type>::data_.size() );
    ArrayBase<type>::add( x , n );
    last_.push_back( ArrayBase<type>::data_.size() );
  }

  void remove( index_t k0 )
  {
    index_t kshift = nv(k0);
    for (index_t k=k0+1;k<nb();k++)
    {
      first_[k] -= kshift;
      last_[k]  -= kshift;
    }
    ArrayBase<type>::remove( first_[k0] , last_[k0] );
    first_.erase( first_.begin() + k0 );
    last_.erase( last_.begin() + k0 );
  }

  std::vector<type> get( index_t k ) const
  {
    std::vector<type> elem;
    for (index_t j=first_[k];j<last_[k];j++)
      elem.push_back(ArrayBase<type>::data_[j]);
    return elem;
  }

  index_t nv( const index_t k ) const
  {
    assert( k < nb() );
    return last_[k] -first_[k];
  }

  index_t nb() const { return first_.size(); }

  type& operator()       (index_t k , index_t j )
    { return ArrayBase<type>::data_[first_[k]+j]; }
  const type& operator() (index_t k , index_t j ) const
    { return ArrayBase<type>::data_[first_[k]+j]; }

  type* operator() (index_t k)
    { return &ArrayBase<type>::data_[first_[k]]; }

  const type* operator() (index_t k) const
    { return &ArrayBase<type>::data_[first_[k]]; }

  type* operator[]       (index_t k)
    { return &ArrayBase<type>::data_[first_[k]]; }

  const type* operator[] (index_t k) const
    { return &ArrayBase<type>::data_[first_[k]]; }

  void clear()
  {
    first_.clear();
    last_.clear();
    ArrayBase<type>::clear();
  }

private:
  std::vector<index_t> first_;
  std::vector<index_t> last_;
};

template<typename type>
class Array<Rectangular,type> : public ArrayBase<type>
{
public:
  Array( coord_t rank ) : rank_(rank) {}

  void set_rank( coord_t rank ) { rank_ = rank; }
  coord_t rank() const { return rank_; }

  void add( const type* x ) { ArrayBase<type>::add(x,rank_); }
  void add( const type* x , index_t n ) { ArrayBase<type>::add(x,rank_); }

  void remove( index_t k0 )
  {
    ArrayBase<type>::remove( k0*rank_ , (k0+1)*rank_ );
  }

  type& operator()       (index_t k , index_t j )
    { return ArrayBase<type>::data_[k*rank_+j]; }
  const type& operator() (index_t k , index_t j ) const
    { return ArrayBase<type>::data_[k*rank_+j]; }

  type* operator[]       (index_t k)
    { return &ArrayBase<type>::data_[k*rank_]; }

  const type* operator[] (index_t k) const
    { return &ArrayBase<type>::data_[k*rank_]; }

  void set( index_t k , const type* x )
  {
    for (coord_t j=0;j<rank_;j++)
      ArrayBase<type>::data_[k*rank_+j] = x[j];
  }

  index_t nv( const index_t ) const
  {
    return rank_;
  }

  index_t nb() const { return ArrayBase<type>::data_.size()/rank_; }

  void clear()
  {
    ArrayBase<type>::clear();
    rank_ = 0;
  }

private:
  coord_t rank_;
};

template<typename type>
class Array<Simple,type> : public ArrayBase<type>
{
public:
  index_t nb() const { return ArrayBase<type>::data_.size(); }

  void add( const type& x ) { ArrayBase<type>::add(&x,1); }

  void remove( index_t k0 )
    { ArrayBase<type>::remove( k0, k0+1 ); }

  void insert( index_t loc , const type& value )
    { ArrayBase<type>::insert(loc,&value,1); }

  type operator[] (index_t k) const
    { return ArrayBase<type>::data_[k]; }

  type&       operator[] (index_t k)
    { return ArrayBase<type>::data_[k]; }

  void set( index_t k , const type& x )
    { ArrayBase<type>::data_[k] = x; }

  void clear()
  {
    ArrayBase<type>::clear();
  }
};

#endif

} // luna

#endif
