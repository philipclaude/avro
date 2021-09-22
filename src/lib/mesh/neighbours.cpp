//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "mesh/neighbours.h"
#include "mesh/topology.h"

namespace avro
{

template<typename type>
index_t
Cache<type>::element( const index_t position ) const
{
  if (position>=element_.size()) print();
  avro_assert_msg( position<element_.size() , "pos = %lu, |cache| = %lu",
                    position,element_.size() );
  return element_[position];
}

template<typename type>
index_t
Cache<type>::index( const index_t position ) const
{
  avro_assert_msg( position<index_.size() , "pos = %lu, |cache| = %lu",
                    position,index_.size() );
  return index_[position];
}

template<typename type>
typename Cache<type>::const_iter
Cache<type>::contains( const std::string& s ) const
{
  return hash_.find(s);
}

template<typename type>
void
Cache<type>::insert( const std::string& s , const index_t elem , const index_t j )
{
  hash_[s] = element_.size();
  element_.push_back(elem);
  index_.push_back(j);
}

template<typename type>
void
Cache<type>::remove( const const_iter& f )
{
  index_t k = f->second;
  hash_.erase( f );
  element_.erase( element_.begin() + k );
  index_.erase( index_.begin() + k );
  std::unordered_map<std::string,index_t>::iterator it;
  for (it=hash_.begin();it!=hash_.end();it++)
  {
    if (it->second>k)
    {
      it->second--;
    }
  }
}

template<typename type>
void
Cache<type>::shift( const index_t elem , const int shift )
{
  for (index_t k=0;k<nb();k++)
  {
    if (element_[k]>elem)
      element_[k] += shift;
  }
}

template<typename type>
std::string
Cache<type>::generate( const std::vector<index_t>& f ) const
{
  if (f.size()==1) return std::to_string(f[0]);
  if (f.size()==2) return std::to_string(f[0])+"|"+std::to_string(f[1]);
  if (f.size()==3) return std::to_string(f[0])+"|"+std::to_string(f[1])+"|"+std::to_string(f[2]);
  if (f.size()==4) return std::to_string(f[0])+"|"+std::to_string(f[1])+"|"+std::to_string(f[2])+"|"+std::to_string(f[3]);
  print_inline(f,"unknown facet dimension");
  avro_implement;
  return "";
}

template<typename type>
std::string
Cache<type>::generate( const index_t element, const index_t j )
{
  // use the nodes to identify the facets (not high-order nodes for that type)
  topology_.facet(element,j,buffer_);
  std::sort( buffer_.begin() , buffer_.end() );
  return generate(buffer_);
}

template<typename type>
template<typename type2>
void
Cache<type>::copy( Cache<type2>& destination )
{
  std::vector<index_t>& e = destination.element();
  e.assign( element_.begin() , element_.end() );

  std::vector<index_t>& i = destination.index();
  i.assign( index_.begin() , index_.end() );

  std::unordered_map<std::string,index_t>& h = destination.hashmap();
  h = hash_;
}

template<typename type>
void
Cache<type>::clear()
{
  hash_.clear();
  element_.clear();
  index_.clear();
}

template<typename type>
void
Cache<type>::print() const
{
  std::unordered_map<std::string,index_t>::const_iterator it;
  for (it=hash_.begin();it!=hash_.end();it++)
  {
    printf("facet[%s] at position %lu, elem = %lu, idx = %lu\n",
            it->first.c_str(),it->second,element(it->second),index(it->second));
  }
}

template<typename type>
Neighbours<type>::Neighbours( Topology<type>& _topology ) :
  topology_(_topology),
  nfacets_(topology_.number()+1),
  facet_( topology_.number() , 0 ),
  cache_(topology_),
  computed_(false),
  fromscratch_(false)
{
}

template<>
void
Neighbours<Simplex>::compute( index_t* elements0 , index_t nelements )
{

  if (computed_)
  {
    printf("neighbours are already computed.\n");
    return;
  }
  cache_.clear();

  std::vector<index_t> elements;
  if (elements0==NULL)
  {
    nelements = topology_.nb();
    elements  = linspace( topology_.nb() );
  }
  else
  {
    elements.resize( nelements );
    for (index_t k=0;k<nelements;k++)
      elements[k] = elements0[k];
  }

  // even though we may be computing a subset of the neighbours
  // they need to be in the correct element location for later retrieval
  neighbours_.resize( topology_.nb()*nfacets_ );
  std::fill( neighbours_.begin() , neighbours_.end() , -1 );

  // add the facets of each element
  for (index_t k=0;k<nelements;k++)
    addElement( elements[k] );

  computed_ = true;

  if (fromscratch_)
    cache_.clear();
}

template<>
void
Neighbours<Polytope>::compute( index_t* elements0 , index_t nelements )
{
  printf("convex polytope neighbours not currently supported.\n");
  avro_implement;
}

template<typename type>
void
Neighbours<type>::findNeighbour( index_t elem0 , index_t j0 )
{

  // generate a string representing this facet
  std::string s = cache_.generate(elem0,j0);

  // search the cache for this facet
  //int position = cache_.contains(s);
  typename Cache<type>::const_iter position = cache_.contains(s);

  if (position==cache_.hashmap().end())
  {
    // cache does not contain, add
    cache_.insert( s , elem0 , j0 );

    //printf("inserting facet %s\n",s.c_str());
  }
  else
  {
    index_t elem1 = cache_.element( position->second );
    index_t j1 = cache_.index( position->second );

    // cache contains this facet, set the neighbours
    neighbours_[ nfacets_*elem0 + j0 ] = elem1;
    neighbours_[ nfacets_*elem1 + j1 ] = elem0;

    // remove the facet from the cache
    if (!fromscratch_)
      cache_.remove( position );

  }
}

template<typename type>
void
Neighbours<type>::enlarge( const index_t nelem )
{
  for (index_t k=0;k<nelem;k++)
  for (index_t j=0;j<nfacets_;j++)
    neighbours_.push_back(-1);
}

template<typename type>
void
Neighbours<type>::remove( const index_t elem0 , bool erase )
{
  // check if the cache contains the element's facet
  for (index_t j0=0;j0<nfacets_;j0++)
  {
    // get the current neighbours of this facet
    int elem1 = operator()(elem0,j0);
    int j1 = indexofme(elem0,j0);

    // set the opposite neighbour to -1
    if (elem1>=0)
    {
      avro_assert( j1>=0 );
      neighbours_[elem1*nfacets_ +j1] = -1;
    }

    // generate a string representing this facet
    std::string s = cache_.generate(elem0,j0);

    //printf("trying to remove facet %s\n",s.c_str());

    // search the cache for this facet
    //int position = cache_.contains(s);
    typename Cache<type>::const_iter position = cache_.contains(s);

    if (position==cache_.hashmap().end() && elem1>=0)
    {
      // cache does not contain this facet so we add it
      // but this time we add the neighbour data since that is
      // what will be used later (i.e. this element disappears)
      // note that if the opposite element is null (elem1<0) then we don't
      // insert
      cache_.insert( s , elem1 , j1 );
    }
    else if (position!=cache_.hashmap().end())
    {
      // cache contains this facet so we remove it entirely
      cache_.remove( position );
    }

    neighbours_[elem0*nfacets_+j0] = -1;

  } // loop over facets

  if (erase)
  {
    // delete the neighbour data of this element
    neighbours_.erase( neighbours_.begin() +elem0*nfacets_ ,
                       neighbours_.begin() +(elem0+1)*nfacets_ );

    // any element neighbour in the neighbour data greater than elem0
    // gets decremented
    for (index_t k=0;k<neighbours_.size();k++)
    {
      if (neighbours_[k]>int(elem0))
        neighbours_[k]--;
    }

    // anything indexing an element larger than elem0 needs to shift down by 1
    cache_.shift( elem0 , -1 );
  }
  else
  {
    std::fill( neighbours_.begin() +elem0*nfacets_ , neighbours_.begin()+(elem0+1)*nfacets_ , -1 );
  }

}

template<typename type>
Cache<type>&
Neighbours<type>::cache()
{
  return cache_;
}

template<typename type>
template<typename type2>
void
Neighbours<type>::copy( Neighbours<type2>& destination )
{
  std::vector<index_t>& f = destination.facet();
  f.assign( facet_.begin() , facet_.end() );

  std::vector<int>& n = destination.neighbours();
  n.assign( neighbours_.begin() , neighbours_.end() );

  Cache<type2>& c = destination.cache();
  cache_.copy(c);

  destination.computed() = true;
}

template<typename type>
bool
Neighbours<type>::cacheMatches( Topology<type>& bnd ) const
{
  if (bnd.nb()!=cache_.nb()) return false;
  for (index_t k=0;k<bnd.nb();k++)
  {
    std::string s = cache_.generate( bnd.get(k) );
    typename Cache<type>::const_iter position = cache_.contains(s);
    if (position==cache_.hashmap().end()) return false;
  }
  return true;
}

template<typename type>
void
Neighbours<type>::print( const index_t k ) const
{
  printf("element %lu: [",k);
  for (index_t j=0;j<nfacets_;j++)
  {
    //printf(" %d (%d)",neighbours_[k*nfacets_+j],indexofme(k,j));
    printf(" %d",neighbours_[k*nfacets_+j]);
  }
  printf(" ]\n");
}

template<typename type>
void
Neighbours<type>::print() const
{
  for (index_t k=0;k<nb();k++)
  {
    printf("element %lu: [",k);
    for (index_t j=0;j<nfacets_;j++)
    {
      printf(" %d (%d)",neighbours_[k*nfacets_+j],indexofme(k,j));
    }
    printf(" ]\n");
  }
  printf("leftovers:\n");
  cache_.print();
}

template class Neighbours<Simplex>;
template class Neighbours<Polytope>;

template void Neighbours<Simplex>::copy( Neighbours<Simplex>& );

} // avro
