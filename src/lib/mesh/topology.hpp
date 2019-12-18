#include "common/tools.h"

#include "mesh/topology.h"
#include "mesh/points.h"

#include <set>

namespace luna
{

template<typename type>
Topology<type>::Topology( Points& points , coord_t number ) :
  Topology(points,number,1)
{}

template<typename type>
void
Topology<type>::get_edges( std::vector<index_t>& edges ) const
{
  std::vector<index_t> ek;

  std::set< std::pair<index_t,index_t> > table;

  for (index_t k=0;k<nb();k++)
  {
    if (ghost(k)) continue;

    const index_t* v0 = operator()(k);

    // get the edges of this cell
    master_.get_edges( v0 , nv(k) , ek );

    // add the edges
    for (index_t j=0;j<ek.size()/2;j++)
    {
      index_t p0 = ek[2*j];
      index_t p1 = ek[2*j+1];

      if (p0<points_.nb_ghost() || p1<points_.nb_ghost())
        continue;

      if (p0>p1) std::swap(p0,p1);
      std::pair<index_t,index_t> E = std::pair<index_t,index_t>(p0,p1);
      if (table.find(E)==table.end())
      {
        table.insert(E);
        edges.push_back(p0);
        edges.push_back(p1);
      }
    }
  }
}

template<typename type>
void
Topology<type>::get_triangles( std::vector<index_t>& triangles ) const
{
  luna_assert( number_==2 );
  triangles.clear();
  for (index_t k=0;k<nb();k++)
  {
    if (ghost(k)) continue;

    for (index_t j=0;j<this->nv(k);j++)
      triangles.push_back( this->operator()(k,j) );
  }
}

template<typename type>
void
Topology<type>::facet( const index_t k , const index_t j ,
                       std::vector<index_t>& f ) const
{
  master_.facet( operator()(k) , j , f );
}

template<typename type>
void
Topology<type>::get_boundary( Topology<type>& boundary ) const
{
  luna_implement;
}

template<typename type>
void
Topology<type>::get_elements( Topology<type>& topology ) const
{
  for (index_t k=0;k<nb_children();k++)
  {
    if (this->child(k).number()!=topology.number())
      continue;
    this->child(k).get_elements(topology);
  }

  if (topology.number()==number())
  {
    for (index_t k=0;k<nb();k++)
      topology.add( (*this)(k) , nv(k) );
  }
}

template<typename type>
bool
Topology<type>::has( index_t k , index_t value ) const
{
  const index_t *elem = (*this)(k);
  for (index_t j=0;j<nv(k);j++)
  {
    if (elem[j]==value)
      return true;
  }
  return false;
}

template<typename type>
void
Topology<type>::all_with( const std::vector<index_t>& f , std::vector<index_t>& elems ) const
{
  elems.clear();
  for (index_t k=0;k<nb();k++)
  {
    // loop through s to see if this element contains it
    bool ok = false;
    for (index_t j=0;j<f.size();j++)
    {
      ok = false;
      for (index_t i=0;i<nv(k);i++)
      {
        if ((*this)(k,i)==f[j])
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
Topology<type>::get_elem( index_t k , std::vector<const real_t*>& X ) const
{
  X.resize( nv(k) );
  for (index_t j=0;j<nv(k);j++)
    X[j] = points_[(*this)(k,j)];
}

template<typename type>
void
Topology<type>::get_elem( index_t k , std::vector<real_t*>& X ) const
{
  X.resize( nv(k) );
  for (index_t j=0;j<nv(k);j++)
    X[j] = points_[(*this)(k,j)];
}

template<typename type>
void
Topology<type>::orient( real_t* q )
{
  if (q!=NULL)
    luna_assert( number_+1 == points_.dim() );

  for (index_t k=0;k<nb();k++)
    orient( operator()(k) , nv(k) , q );
}

template<typename type>
real_t
Topology<type>::volume() const
{
  real_t v = 0.;
  for (index_t k=0;k<nb();k++)
  {
    v += master_.volume( points_ , operator()(k) , nv(k) );
  }
  return v;
}

template<typename type>
void
Topology<type>::remove_point( const index_t k )
{
  points_.remove(k);

  // decrement any indices higher than the original index
  for (index_t i=0;i<data_.size();i++)
  {
    if (data_[i]>k)
      data_[i]--;
  }
}

} // luna
