#include "common/array.h"
#include "common/tools.h"

#include "mesh/facets.h"
#include "mesh/topology.h"

#include <algorithm>

namespace luma
{

Facets::Facets( const Topology<Simplex>& topology ) :
  topology_(topology),
  elem2facet_(TableLayout_Rectangular)
{}

index_t
Facets::opposite( const index_t f , short side ) const
{
  if (side==0) return topology_( side0_[f] , indx0_[f] );
  if (side==1) return topology_( side1_[f] , indx1_[f] );
  luma_assert_not_reached;
  return 0;
}

index_t
Facets::opposite( const index_t f ) const
{
  luma_assert( boundary(f) );
  if (side0_[f]<0) return topology_( side1_[f] , indx1_[f] );
  else return topology_(side0_[f] , indx0_[f] );
}

index_t
Facets::opposite( std::vector<index_t>& F ) const
{
  std::sort(F.begin(),F.end());
  std::vector<index_t> f( topology_.number() );
  for (index_t k=0;k<nb();k++)
  {
    if (!boundary(k)) continue;

    retrieve(k,f);
    std::sort( f.begin() , f.end() );

    // check if it's equal
    bool same = true;
    for (index_t i=0;i<f.size();i++)
    {
      if (f[i]!=F[i])
      {
        same = false;
        break;
      }
    }
    if (!same) continue;

    return opposite( k );
  }
  luma_assert_not_reached;
  return 0;
}

void
Facets::compute( index_t* elements0 , index_t nelem )
{
  clear();

  std::vector<index_t> elements;
  if (elements0==NULL)
  {
    elements = linspace(topology_.nb());
    nelem = topology_.nb();
  }
  else
  {
    elements.resize( nelem );
    for (index_t k=0;k<nelem;k++)
      elements[k] = elements0[k];
  }

  coord_t nf = topology_.number()+1;
  std::vector<index_t> f( topology_.number() );
  elem2facet_.set_rank( nf );

  // even in case a subset of the elements are used, the elem2facet info
  // is sized for all the elements because the indices need to correspond
  // to the correct elements
  elem2facet_.allocate( topology_.nb() );

  // loop through all elements and construct the facets of each
  for (index_t elem=0;elem<nelem;elem++)
  {

    index_t k = elements[elem];

    // loop through facets
    for (coord_t j=0;j<nf;j++)
    {
      facet( k , j , f );
      HashableElement s(f.data(),f.size());
      int id;
      id = __contains__(s);
      if (id<0)
      {
        id = nb();

        // add the facet and the relevant information
        facet_.push_back(s);
        _facet_.insert( std::make_pair(s,_facet_.size()) );
        side0_.push_back(k);
        indx0_.push_back(j);

        side1_.push_back(-1); // out of bounds
        indx1_.push_back(-1); // out of bounds

        // TODO cache the facet
      }
      else
      {
        // set the relevant information
        side1_[id] = k;
        indx1_[id] = j;
      }

      // save the element-to-facet information
      elem2facet_(k,j) = id;
    }
  }
}

void
Facets::facet( const index_t* v , const index_t j , std::vector<index_t>& f ) const
{
  index_t i = 0;
  for (coord_t d=0;d<f.size()+1;d++)
  {
    if (d==j) continue;
    f[i++] = v[d];
  }
}

index_t
Facets::elem2facet( const index_t k , const index_t j ) const
{
  return elem2facet_(k,j);
}

void
Facets::facet( index_t k , const index_t j , std::vector<index_t>& f ) const
{
  facet( topology_(k) , j , f );
}

int
Facets::__contains__( HashableElement& s ) const
{
  std::unordered_map<HashableElement,index_t>::const_iterator it = _facet_.find( s );
  if (it==_facet_.end()) return -1;
  return it->second;
}

index_t
Facets::neighbour( const index_t k , const index_t j ) const
{
  const index_t f = elem2facet_(k,j);
  if ( side0_[f] == int(k) ) return side1_[f];
  luma_assert( side1_[f]==int(k) );
  return side0_[f];
}

bool
Facets::check() const
{
  const index_t nf = topology_.number()+1;

  bool result = true;

  for (index_t k=0;k<nb();k++)
  {
    if (side0_[k]<0)
    {
      printf("side0_[%lu] = %d\n",k,side0_[k]);
      result = false;
    }
    if (indx0_[k]<0)
    {
      printf("indx0_[%lu] = %d\n",k,indx0_[k]);
      result = false;
    }
    if (side1_[k]<0)
    {
      printf("side1_[%lu] = %d\n",k,side1_[k]);
      result = false;
    }
    if (indx1_[k]<0)
    {
      printf("indx1_[%lu] = %d\n",k,indx1_[k]);
      result = false;
    }
  }

  luma_assert( elem2facet_.nb()==topology_.nb() );
  for (index_t k=0;k<topology_.nb();k++)
  {
    for (index_t j=0;j<nf;j++)
    {
      if (elem2facet_(k,j)<0)
        result = false;
    }
  }

  return result;
}

int
Facets::incache( const HashableElement& s ) const
{
  for (index_t k=0;k<cache_.size();k++)
  {
    if (facet_[cache_[k]]==s)
      return k;
  }
  return -1;
}

void
Facets::retrieve( const index_t fid , std::vector<index_t>& f ) const
{
  // do NOT decompose the string label
  // look at either side0 or side1 and get the elements facet
  // reason: one of them might be a boundary facet
  if (side0_[fid]>=0)
  {
    facet( side0_[fid] , indx0_[fid] ,  f );
    return;
  }
  luma_assert_msg( side1_[fid]>=0 , "this facet (%lu) has no neighbours!" , fid );

  facet( side1_[fid] , indx1_[fid] , f );
}

void
Facets::retrieve( const index_t fid , std::vector<int>& T ) const
{
  // do NOT decompose the string label
  // look at either side0 or side1 and get the elements facet
  // reason: one of them might be a boundary facet
  std::vector<index_t> f( T.size()-1 );
  if (side0_[fid]>=0)
  {
    // fill T with the side0 element but replace indx0 with -1
    for (index_t j=0;j<T.size();j++)
      T[j] = topology_( side0_[fid] , j );
    T[ indx0_[fid] ] = -1;
    return;
  }
  luma_assert_msg( side1_[fid]>=0 , "this facet (%lu) has no neighbours!" , fid );

  // fill T with the side0 element but replace indx0 with -1
  for (index_t j=0;j<T.size();j++)
    T[j] = topology_( side0_[fid] , j );
  T[ indx0_[fid] ] = -1;
}


void
Facets::clear()
{
  side0_.clear();
  indx0_.clear();
  side1_.clear();
  indx1_.clear();
  elem2facet_.clear();
  facet_.clear();
  _facet_.clear();
}

void
Facets::print() const
{
  for (index_t k=0;k<nb();k++)
  {
    printf("facet[%lu]: side0 = %d (%d), side1 = %d (%d)\n",
      k,side0_[k],indx0_[k],side1_[k],indx1_[k]);
  }
  elem2facet_.print();
}

} // luma
