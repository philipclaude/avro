#ifndef LUNA_LIB_MESH_NEIGHBOURS_H_
#define LUNA_LIB_MESH_NEIGHBOURS_H_

#include "common/types.h"
#include "common/tools.h"

#include <algorithm>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace luna
{

template<typename type> class Topology;

// implements a cache for the facets which only point to one neighbour
template<typename type>
class Cache
{
public:
  Cache( Topology<type>& _topology ) :
    topology_(_topology),
    buffer_(topology_.number(),0)
 {}

   typedef std::unordered_map<std::string,index_t>::const_iterator const_iter;

  index_t nb() const { return element_.size(); }
  index_t element( const index_t position ) const;
  index_t index( const index_t position ) const;

  const_iter contains( const std::string& s ) const;
  void insert( const std::string& s , const index_t elem , const index_t j );
  void remove( const const_iter& it );
  void shift( const index_t elem , const int shift );

  std::string generate( const std::vector<index_t>& f ) const;
  std::string generate( const index_t element, const index_t j );

  template<typename type2>
  void copy( Cache<type2>& destination );

  void clear();

  std::vector<index_t>& element() { return element_; }
  std::vector<index_t>& index() { return index_; }
  std::unordered_map<std::string,index_t>& hashmap() { return hash_; }
  const std::unordered_map<std::string,index_t>& hashmap() const { return hash_; }

  void print() const;

private:
  Topology<type>& topology_;

  std::vector<index_t> element_;
  std::vector<index_t> index_;

  std::unordered_map<std::string,index_t> hash_;
  std::vector<index_t> buffer_;


};

// templated by element type to traverse through Simplex and CurvilinearSimplex
// storage only requires number+1 neighbours per element in either of
// these cases
template<typename type>
class Neighbours
{
private:


public:
  Neighbours( Topology<type>& _topology );

  void allocate();
  void compute( index_t* elements0=NULL , index_t nelements=0 );

  void findNeighbour( index_t elem , index_t j );

  void addElement( index_t elem )
  {
    for (index_t j=0;j<nfacets_;j++)
      findNeighbour(elem,j);
  }

  void enlarge( const index_t nelem );

  void facet( const index_t k , const index_t j , std::vector<index_t>& f )
  {
    index_t count = 0;
    f.resize(nfacets_-1);
    for (index_t i=0;i<nfacets_;i++)
    {
      if (i==j) continue;
      f[count++] = topology_(k,i);
    }
  }

  index_t nfacets() const { return nfacets_; }

  void remove( const index_t k , bool erase=true );

  int indexofme( const index_t k , const index_t j ) const
  {
    // return the index of element k in its neighbour stored at j
    int n = operator()(k,j);
    if (n<0) return n;
    for (index_t i=0;i<nfacets_;i++)
      if (operator()(index_t(n),i)==int(k))
        return i;
    luna_assert_not_reached;
    return 0;
  }

  int operator() ( const index_t k , const index_t j ) const
    { return neighbours_[ k*nfacets_+j ]; }

  int opposite( const index_t k0 , const index_t k1 ) const
  {
    for (index_t j=0;j<nfacets_;j++)
    {
      if (operator()(k0,j)==int(k1))
        return int(topology_(k0,j));
    }
    return -1;
  }

  short oppositeIndex( const index_t k0 , const index_t k1 ) const
  {
    for (coord_t j=0;j<nfacets_;j++)
    {
      if (operator()(k0,j)==int(k1))
        return j;
    }
    luna_assert_not_reached;
    return -1;
  }

  index_t nb() const { return neighbours_.size()/nfacets_; }

  void print() const;
  void print( const index_t k ) const;

  Cache<type>& cache();

  bool cacheMatches( Topology<type>& bnd ) const;

  template<typename type2>
  void copy( Neighbours<type2>& destination );

  void forceCompute() { computed_ = false; }
  bool& computed() { return computed_; }

  bool& fromscratch() { return fromscratch_; }

  std::vector<index_t>& facet() { return facet_; }
  std::vector<int>& neighbours() { return neighbours_; }

private:
  Topology<type>& topology_;
  const coord_t nfacets_;
  std::vector<index_t> facet_;

  Cache<type> cache_;
  std::vector<int> neighbours_;

  bool computed_;

  bool fromscratch_;
};

} // luna

#endif
