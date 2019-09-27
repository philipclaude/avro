#ifndef URSA_COMMON_DATA_H_
#define URSA_COMMON_DATA_H_

#include "common/types.h"

#include <map>
#include <string>
#include <vector>

namespace ursa
{

template<typename type>
class Data
{

public:
  Data() : sorted_(false), preallocated_(false) {}
  Data( const bool _sorted ) : sorted_(_sorted) {}

  // data & information functions
  std::vector<type> elements() const { return elements_; }
  std::vector<index_t> first() const { return first_; }
  std::vector<index_t> last() const { return last_; }
  bool sorted() const { return sorted_; }
  bool offset() const { return offset_; }
  void setOffset( const bool x ) { offset_ = x; }
  void setSorted( const bool x ) { sorted_ = x; }

  void sort();

  // allocators
  void allocate( const index_t n , const index_t size );
  void allocate( const index_t n , const std::vector<index_t> sizes );

  // data addition & modification
  void add( std::vector<type> elem );
  void add( type* d0 , const index_t nd );
  void add( const type* d0 , const index_t nd ) { __add__(d0,nd); }
  void addto( const index_t k , const index_t value );
  void offsetBy( const type offset );
  void incrementIfGreater( const type offset , const type d0 );
  void decrementIfGreater( const type offset , const type d0 );
  void decrement( const type offset );
  void set( const index_t k , const index_t j , const index_t value );
  void setall( const type& value );

  void mapData( std::map<type,type>& dmap );

  // set functions
  bool has( const index_t k , const type value ) const;
  bool has( const type value ) const;
  index_t countOccurrencesOf( const type n );
  index_t elementsWith( const type n , std::vector<index_t>& elems ) const;
  void allWithSubset( const std::vector<type>& sub , std::vector<index_t>& elems ) const;
  void dataOfElems( const std::vector<index_t>& elems , std::vector<index_t>& data );
  void closure( const std::vector<index_t>& a0 , std::vector<index_t>& a , std::vector<index_t>& N , Data<index_t>& da );
  void remove( const index_t k0 );
  void remove( const index_t k0 , const index_t j );
  void replace( const index_t k0 , type* v1 , const index_t nv1 );
  void operate( const std::vector<index_t>& subtractions , Data<type>& additions );
  bool contains( type* d , index_t nd ) const;
  index_t cardinality( type*d , index_t nd ) const;

  // size information
  index_t nb() const { return first_.size(); }
  index_t nv( const index_t k ) const { return index_t(last_[k]-first_[k]); }
  type max() const;

  // retrieval
  type* operator() ( const index_t k );
  type operator() ( const index_t k , const index_t j ) const;
  type& operator() ( const index_t k , const index_t j );
  std::vector<type> get( const index_t k ) const;

  // utilities
  void clear();
  void printData( const index_t nt=0 , const std::string& name=std::string() , const std::string& prefix=std::string() ) const;

  void reserve( index_t nelem , index_t ndof_per_elem );

  bool vacant( const index_t k ) const { return vacant_[k]; }
  void preallocate( index_t nelem , index_t nv_per_elem );
  int findvacant();

private:

  void __add__( type* d0 , const index_t nd )
  {
    if (sorted_) std::sort( d0 , d0+nd );
    first_.push_back( elements_.size() );
    for (index_t j=0;j<nd;j++)
      elements_.push_back( d0[j] );
    last_.push_back( elements_.size() );
  }

  void __add__( const type* d0 , const index_t nd )
  {
    // no sorting
    first_.push_back( elements_.size() );
    for (index_t j=0;j<nd;j++)
      elements_.push_back( d0[j] );
    last_.push_back( elements_.size() );
  }

protected:
  std::vector<type> elements_;
  std::vector<index_t> first_;
  std::vector<index_t> last_;
  std::vector<bool> vacant_;
  bool offset_;
  bool sorted_;
  bool preallocated_;

};

} // ursa

#endif
