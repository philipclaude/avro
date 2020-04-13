#ifndef avro_LIB_TABLE_H_
#define avro_LIB_TABLE_H_

#include "common/error.h"
#include "common/tools.h"
#include "common/types.h"

#include <algorithm>
#include <vector>

namespace avro
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
  Table( TableLayoutCategory layout = TableLayout_None ) :
    layout_(layout),
    sorted_(false),
    rank_(0)
  {}

  Table( const TableLayoutCategory& layout , index_t rank ) :
    layout_(layout),
    sorted_(false),
    rank_(rank)
  {}

  bool undefined() const { return layout_==TableLayout_None; }

  void allocate( index_t n )
  {
    avro_assert( layout_ == TableLayout_Rectangular );
    data_.resize( n*rank_ , type(0) );
  }

  void copy( const Table<type>& table )
  {
    layout_ = table.layout();
    sorted_ = table.sorted();
    rank_   = table.rank();

    data_  = table.data();
    first_ = table.first();
    last_  = table.last();
  }

  void set_layout( TableLayoutCategory layout )
    { layout_ = layout; }
  TableLayoutCategory layout() const { return layout_; }

  void set_sorted( bool sorted ) { sorted_ = sorted; }
  bool sorted() const { return sorted_; }

  void set_rank( index_t rank )
  {
    avro_assert( layout_ == TableLayout_Rectangular );
    rank_ = rank;
  }
  index_t rank() const
  {
    //avro_assert( layout_==TableLayout_Rectangular );
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
  {
    avro_assert_msg( location(k)+j < data_.size() , "k = %lu, j = %lu, nb = %lu, nv = %lu",k,j,nb(),nv(k) );
    return data_[location(k)+j];
  }
  const type operator() (index_t k, index_t j) const
  {
      avro_assert_msg( location(k)+j < data_.size() , "k = %lu, j = %lu, nb = %lu, nv = %lu",k,j,nb(),nv(k) );
      return data_[location(k)+j];
  }

  index_t nv( index_t k ) const
  {
    if (layout_==TableLayout_Rectangular) return rank_;
    else if (layout_==TableLayout_Jagged) return last_[k]-first_[k];
    else avro_assert_not_reached;
  }

  index_t nb() const
  {
    if (layout_==TableLayout_Rectangular) return data_.size()/rank_;
    else if (layout_==TableLayout_Jagged) return first_.size();
    else avro_assert_not_reached;
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
    if (layout_==TableLayout_Rectangular)
      avro_assert( n == rank_ );
    if (layout_==TableLayout_Jagged)
      first_.push_back( data_.size() );
    for (index_t j=0;j<n;j++)
      data_.push_back(x[j]);
    if (layout_==TableLayout_Jagged)
      last_.push_back( data_.size() );
  }

  void add( const type* x , index_t n )
  {
    if (layout_==TableLayout_Rectangular)
      avro_assert_msg( n == rank_ , "n = %lu, rank = %lu" , n , rank_ );
    if (layout_==TableLayout_Jagged)
      first_.push_back( data_.size() );
    for (index_t j=0;j<n;j++)
      data_.push_back(x[j]);
    if (layout_==TableLayout_Jagged)
      last_.push_back( data_.size() );
  }

  bool has( index_t k , type x ) const
  {
    for (index_t j=location(k);j<location(k)+nv(k);j++)
      if (data_[j]==x)
        return true;
    return false;
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
    if (layout_==TableLayout_Rectangular)
    {
      remove( k0*rank_ , (k0+1)*rank_ );
    }
    if (layout_==TableLayout_Jagged)
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

  std::vector<index_t> first() const { return first_; }
  std::vector<index_t> last() const { return last_; }

  void remove( index_t start , index_t end )
  {
    data_.erase( data_.begin() + start , data_.begin() + end );
  }

  index_t location( index_t k ) const
  {
    if (layout_==TableLayout_Rectangular) return k*rank_;
    else if (layout_==TableLayout_Jagged) return first_[k];
    else avro_assert_not_reached;
  }

protected:
  TableLayoutCategory layout_;
  bool sorted_;
  index_t rank_;

  std::vector<type>    data_;
  std::vector<index_t> first_;
  std::vector<index_t> last_;

};

} // avro

#endif
