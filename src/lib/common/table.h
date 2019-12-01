#ifndef LUNA_LIB_TABLE_H_
#define LUNA_LIB_TABLE_H_

#include "common/error.h"
#include "common/types.h"

#include <algorithm>
#include <vector>

namespace luna
{

enum TableLayoutCategory
{
  TableLayout_Jagged,
  TableLayout_Rectangular,
  TableLayout_None
};

// look-up table
template<typename type>
class Table
{
public:
  Table( TableLayoutCategory category = TableLayout_None ) :
    category_(category),
    sorted_(false),
    rank_(0)
  {}

  Table( const TableLayoutCategory& category , index_t rank ) :
    category_(category),
    sorted_(false),
    rank_(rank)
  {}

  bool undefined() const { return category_==TableLayout_None; }

  void allocate( index_t n )
  {
    luna_assert( category_ == TableLayout_Rectangular );
    data_.resize( n*rank_ , type(0) );
  }

  void set_category( TableLayoutCategory category )
    { category_ = category; }

  void set_sorted( bool sorted ) { sorted_ = sorted; }

  void set_rank( index_t rank )
  {
    luna_assert( category_ == TableLayout_Rectangular );
    rank_ = rank;
  }

  index_t rank() const
  {
    luna_assert( category_==TableLayout_Rectangular );
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
    if (category_==TableLayout_Rectangular) return rank_;
    else if (category_==TableLayout_Jagged) return last_[k]-first_[k];
    else luna_assert_not_reached;
  }

  index_t nb() const
  {
    if (category_==TableLayout_Rectangular) return data_.size()/rank_;
    else if (category_==TableLayout_Jagged) return first_.size();
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

  void add( type* x , index_t n )
  {
    if (category_==TableLayout_Rectangular)
      luna_assert( n == rank_ );
    if (category_==TableLayout_Jagged)
      first_.push_back( data_.size() );
    for (index_t j=0;j<n;j++)
      data_.push_back(x[j]);
    if (category_==TableLayout_Jagged)
      last_.push_back( data_.size() );
  }

  void add( const type* x , index_t n )
  {
    if (category_==TableLayout_Rectangular)
      luna_assert( n == rank_ );
    if (category_==TableLayout_Jagged)
      first_.push_back( data_.size() );
    for (index_t j=0;j<n;j++)
      data_.push_back(x[j]);
    if (category_==TableLayout_Jagged)
      last_.push_back( data_.size() );
  }

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
    if (category_==TableLayout_Rectangular)
    {
      remove( k0*rank_ , (k0+1)*rank_ );
    }
    if (category_==TableLayout_Jagged)
    {
      index_t kshift = nv(k0);
      for (index_t k=k0+1;k<nb();k++)
      {
        first_[k] -= kshift;
        last_[k]  -= kshift;
      }
      remove( first_[k0] , last_[k0] );
      first_.erase( first_.begin() + k0 );
      last_.erase( last_.begin() + k0 );
    }
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

  void remove( index_t start , index_t end )
  {
    data_.erase( data_.begin() + start , data_.begin() + end );
  }

  index_t location( index_t k ) const
  {
    if (category_==TableLayout_Rectangular) return k*rank_;
    else if (category_==TableLayout_Jagged) return first_[k];
    else luna_assert_not_reached;
  }

  TableLayoutCategory category_;
  bool sorted_;
  index_t rank_;

  std::vector<type>    data_;
  std::vector<index_t> first_;
  std::vector<index_t> last_;

};

} // luna

#endif