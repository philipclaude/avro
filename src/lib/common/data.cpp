#include "common/data.h"
#include "common/tools.h"

#include "numerics/matrix.h"

#include <algorithm>

namespace luna
{

template<typename type>
void
Data<type>::sort()
{
  for (index_t k=0;k<nb();k++)
    std::sort( elements_.begin()+first_[k] , elements_.begin()+last_[k]-1 );
  sorted_ = true;
}

template<typename type>
void
Data<type>::allocate( const index_t n , const index_t size )
{
  for (index_t k=0;k<n;k++)
  {
    first_.push_back( elements_.size() );
    for (index_t j=0;j<size;j++)
      elements_.push_back(type(0));
    last_.push_back( elements_.size() );
  }
}

template<typename type>
void
Data<type>::allocate( const index_t n , const std::vector<index_t> sizes )
{
  for (index_t k=0;k<n;k++)
  {
    first_.push_back( elements_.size() );
    for (index_t j=0;j<sizes[k];j++)
      elements_.push_back(type(0));
    last_.push_back( elements_.size() );
  }
}

template<typename type>
type*
Data<type>::operator() ( const index_t k )
{
  return &elements_[first_[k]];
}

template<typename type>
const type*
Data<type>::operator() ( const index_t k ) const
{
  return &elements_[first_[k]];
}

template<typename type>
type
Data<type>::operator() ( const index_t k , const index_t j ) const
{
  luna_assert_msg( k < first_.size() && k < last_.size() ,
                  "k = %d, |first| = %d, |last| = %d\n",
                  int(k),int(first_[k]),int(last_[k]) );
  luna_assert_msg( first_[k] +j < elements_.size() ,
                   "first_[%d] = %d, j = %d, |elements| = %d\n",
                  int(k),int(first_[k]),int(j),int(elements_.size()) );
  return elements_[ first_[k] +j ];
}

template<typename type>
type&
Data<type>::operator() ( const index_t k , const index_t j )
{
  luna_assert_msg( k < first_.size() && k < last_.size() ,
                  "k = %d, |first| = %d, |last| = %d\n",
                  int(k),int(first_[k]),int(last_[k]) );
  luna_assert_msg( first_[k] +j < elements_.size() ,
                  "first_[%d] = %d, j = %d, |elements| = %d\n",
                  int(k),int(first_[k]),int(j),int(elements_.size()) );
  return elements_[ first_[k] +j ];
}

template<typename type>
std::vector<type>
Data<type>::get( const index_t k ) const
{
  std::vector<type> elem;
  for (index_t j=first_[k];j<last_[k];j++)
    elem.push_back(elements_[j]);
  return elem;
}

template<typename type>
bool
Data<type>::has( const index_t k , const type value ) const
{
  for (index_t j=first_[k];j<last_[k];j++)
  {
    if (elements_[j]==value)
      return true;
  }
  return false;
}

template<typename type>
bool
Data<type>::has( const type value ) const
{
  for (index_t k=0;k<elements_.size();k++)
    if (elements_[k]==value) return true;
  return false;
}

template<typename type>
void
Data<type>::add( std::vector<type> elem )
{
  __add__(elem.data(),elem.size());
}

template<typename type>
void
Data<type>::add( type* d0 , const index_t nd )
{
  __add__(d0,nd);
}

template<typename type>
void
Data<type>::clear()
{
  elements_.clear();
  first_.clear();
  last_.clear();
  offset_ = false;
}

template<>
void
Data<index_t>::allwith( const std::vector<index_t>& sub , std::vector<index_t>& elems ) const
{
  elems.clear();
  for (index_t k=0;k<nb();k++)
  {
    // loop through s to see if this element contains it
    bool ok = false;
    for (index_t j=0;j<sub.size();j++)
    {
      ok = false;
      for (index_t i=0;i<nv(k);i++)
      {
        if (operator()(k,i)==sub[j])
        {
          // the element has this index
          ok = true;
          break;
        }
      }
      // this index was not found so the element does not have this facet
      if (!ok) break;
      //else printf("elem %lu contains index %lu\n",k,sub[j]);
    }

    // if we made it here with ok being true, then the element has the facet
    if (ok) elems.push_back(k);
  }
}

template<typename type>
void
Data<type>::remove( const index_t k0 )
{
  luna_assert_msg( first_.size()==nb() ,
              "nb = %d but first size = %d", int(nb()) , int(first_.size()) );
  luna_assert_msg( last_.size()==nb() ,
              "nb = %d but last size = %d", int(nb()) , int(last_.size()) );
  index_t kshift = nv(k0);
  for (index_t k=k0+1;k<nb();k++)
  {
    first_[k] -= kshift;
    last_[k]  -= kshift;
  }
  elements_.erase( elements_.begin() +first_[k0] , elements_.begin()+last_[k0] );
  first_.erase( first_.begin() +k0 );
  last_.erase( last_.begin() +k0 );
}

template<typename type>
bool
Data<type>::contains( type* d , index_t nd ) const
{
  bool has;
  for (index_t k=0;k<nb();k++)
  {
    has = true;
    if (nd!=nv(k)) return false;
    for (index_t j=0;j<nv(k);j++)
    {
      if (d[j]!=operator()(k,j))
      {
        has = false;
        break;
      }
    }
    if (has) return true;
  }
  return false;
}

template<typename type>
index_t
Data<type>::cardinality( type* d0 , index_t nd ) const
{
  std::vector<type> d(d0,d0+nd);
  std::sort(d.begin(),d.end());
  index_t count = 0;
  bool has;
  std::vector<type> tk(nd);
  for (index_t k=0;k<nb();k++)
  {
    has = true;
    for (index_t j=0;j<nd;j++)
      tk[j] = operator()(k,j);
    if (nd!=tk.size()) return false;
    std::sort(tk.begin(),tk.end());
    for (index_t j=0;j<tk.size();j++)
    {
      if (d[j]!=tk[j])
      {
        has = false;
        break;
      }
    }
    if (has) count++;
  }
  return count;
}

template<typename type>
type
Data<type>::max() const
{
  if (elements_.size()==0) return type(0);
  return *std::max_element( elements_.begin() , elements_.end() );
}

template<>
void
Data<index_t>::map( std::map<index_t,index_t>& dmap )
{
  for (index_t k=0;k<elements_.size();k++)
  {
    elements_[k] = dmap[elements_[k]];
  }
}

template<typename type>
void
Data<type>::reserve( index_t nelem , index_t ndof_per_elem )
{
  first_.reserve( nelem );
  last_.reserve( nelem );
  elements_.reserve( nelem*ndof_per_elem );
}

template<typename type>
void
Data<type>::setall( const type& value )
{
  std::fill( elements_.begin() , elements_.end() , value );
}

template<typename type>
void
Data<type>::print( const index_t nt , const std::string& name , const std::string& prefix ) const
{
  for (index_t k=0;k<nt;k++) printf("  ");
  if (!name.empty()) printf("%s: (type=%s)\n",name.c_str(),typeid(type).name());
  printf("sorted = %s, offset = %s\n",(sorted_)?"true":"false",(offset_)?"true":"false");
  for (index_t k=0;k<nb();k++)
  {
    std::vector<type> d(nv(k));
    for (index_t j=0;j<nv(k);j++) d[j] = operator()(k,j);
    if (prefix.empty()) print_inline( d , "data" , k );
    else print_inline( d, prefix , k );
  }
}

template class Data<index_t>;
template class Data<real_t>;
template class Data<int>;
template class Data<std::vector<real_t>>;

} // luna
