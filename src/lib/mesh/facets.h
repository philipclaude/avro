#ifndef LUNA_LIB_MESH_FACETS_H_
#define LUNA_LIB_MESH_FACETS_H_

#include "common/table.h"
#include "common/types.h"

#include <algorithm>
#include <functional>
#include <string>
#include <vector>
#include <unordered_map>

namespace luna
{

class Simplex;
template <typename type> class Topology;

class HashableElement
{
public:
  HashableElement( const index_t* v , const index_t nv ) :
    elem_(v,v+nv)
  {
    std::sort(elem_.begin(),elem_.end());
  }

  bool operator==( const HashableElement& elem ) const
  {
    bool equal = true;
    for (std::size_t i = 0; (i < elem_.size()) && equal; i++)
      equal = equal && (elem_[i] == elem.elem()[i]);
    return equal;
  }

  bool operator!=( const HashableElement& elem ) const
    { return !operator==(elem); }


  const std::vector<index_t>& elem() const { return elem_; }

private:
  std::vector<index_t> elem_;
};

} // luna

namespace std
{

template<> struct hash<luna::HashableElement>
{
  std::size_t operator() (luna::HashableElement const& elem) const noexcept
  {
    std::size_t seed = 0;
    for(auto& i : elem.elem())
        seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    /*
    boost::hash_combine(seed, elem.topo);
    for (int n : elem.sortedNodes())
      boost::hash_combine(seed, n);
    */
    return seed;
  }
};

} // namespace std

namespace luna
{

  inline std::size_t hash_value(const HashableElement& elem)
  {
    std::size_t seed = 0;
    /*
    boost::hash_combine(seed, elem.topo);
    for (int n : elem.sortedNodes())
      boost::hash_combine(seed, n);
    */
    return seed;
  }

class Facets
{
public:
  Facets( const Topology<Simplex>& topology );

  void compute( index_t* elem0=NULL , index_t nelem=0 );

  index_t side0( const index_t f ) const { return side0_[f]; }
  index_t side1( const index_t f ) const { return side1_[f]; }

  index_t opposite( const index_t f , short side ) const;
  index_t opposite( const index_t f ) const;
  index_t opposite( std::vector<index_t>& F ) const;

  void facet( index_t k , index_t j , std::vector<index_t>& f ) const;
  void facet( const index_t* k , const index_t j , std::vector<index_t>& f ) const;
  index_t neighbour( const index_t k , const index_t j ) const;
  index_t elem2facet( const index_t k , const index_t j ) const;

  index_t nb() const { return side0_.size(); }

  void print() const;

  bool check() const;

  void clear();

  index_t nb_boundary() const;
  bool boundary( const index_t k ) const
  {
    if (side0_[k]<0) return true;
    if (side1_[k]<0) return true;
    return false;
  }

  const HashableElement& facetid( const index_t k ) const { return facet_[k]; }

  void retrieve( const index_t fid , std::vector<index_t>& f ) const;
  void retrieve( const index_t fid , std::vector<int>& f ) const;

private:

  int incache( const HashableElement& s ) const;
  void clearCache() { cache_.clear(); }
  std::vector<index_t> cache_;

  const Topology<Simplex>& topology_;

  Table<int> elem2facet_;

  // first side of facet
  std::vector<int>   side0_;
  std::vector<short> indx0_; // index of this facet relative to the element
                                 // equivalent to f = elem2facet_(side0)[indx0]

  // second side of facet
  std::vector<int>   side1_;
  std::vector<short> indx1_;

  // only used in construction, not tracked
  std::vector<HashableElement> facet_;
  std::unordered_map<HashableElement,index_t> _facet_;
  int __contains__( HashableElement& fs ) const;
};

} // luna

#endif
