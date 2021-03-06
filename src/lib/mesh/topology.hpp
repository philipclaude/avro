//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "common/tools.h"

#include "mesh/topology.h"
#include "mesh/points.h"

#include <algorithm>
#include <list>
#include <set>
#include <unordered_set>

namespace avro
{

template<typename type>
Topology<type>::Topology( Points& points , coord_t number ) :
  Topology(points,number,1)
{}

template<typename type>
template<typename Friend_t>
void
Topology<type>::construct( std::shared_ptr<Topology<Friend_t>>& node , Topology<Friend_t>& root ) const
{
  node = std::make_shared<Topology<Friend_t>>(root.points(),number_,element_.order());
}

template<typename type>
void
get_edges_using_ball( const Topology<type>& topology , std::vector<index_t>& edges )
{
  std::vector<index_t> ball;
  std::vector<index_t> N;

  for (index_t k = 0; k < topology.points().nb(); k++) {

    if (k < topology.points().nb_ghost()) continue;
    ball.clear();
    N.clear();

    index_t p = k;
    topology.inverse().ball(p,ball);

    for (index_t j = 0; j < ball.size(); j++) {
      for (index_t i = 0; i < topology.nv(ball[j]); i++) {
        index_t q = topology(ball[j],i);
        if (q < topology.points().nb_ghost()) continue;
        if (p == q) continue;
        if (p > q) continue;
        N.push_back(q);
      }
    }
    uniquify(N);

    for (index_t j = 0; j < N.size(); j++) {
      edges.push_back(p);
      edges.push_back(N[j]);
    }
  }
}

template<typename type>
void
Topology<type>::get_edges( std::vector<index_t>& edges ) const
{

  if (inverse_.created()) {
    //get_edges_using_ball(*this,edges);
    //return;
  }

  std::vector<index_t> ek;

  std::set< std::pair<index_t,index_t> > table;

  for (index_t k=0;k<nb();k++)
  {
    if (ghost(k)) continue;

    const index_t* v0 = operator()(k);

    // get the edges of this cell
    element_.get_edges( v0 , nv(k) , ek );

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
Topology<type>::get_topologies( std::vector<const TopologyBase*>& c ) const
{
  for (index_t k=0;k<nb_children();k++)
  {
    c.push_back(this->child_ptr(k));
    this->child(k).get_topologies(c);
  }
  uniquify(c);
}

#if 0

template<typename type>
void
Topology<type>::triangulate( Topology<Simplex>& triangulation ) const
{

  if (number_==0)
  {
    // add the simplices
    for (index_t k=0;k<nb();k++)
      triangulation.add( operator()(k,0) );
  }

  // loop through the cells
  for (index_t k=0;k<nb();k++)
  {

    if (ghost(k)) continue;

    // record the number of vertices and simplices we currently have
    index_t np0 = triangulation.points().nb();
    index_t nt0 = triangulation.nb();

    // ask the  to triangulate -- it needs these vertices to compute the
    // geometry of the vertices it creates
    element_.triangulate( triangulation , (*this)(k) , nv(k) );

    // loop through the created simplices
    for (index_t j=nt0;j<triangulation.nb();j++)
    {
      avro_assert( triangulation.nv(j)==index_t(triangulation.number()+1) );
      for (index_t i=0;i<triangulation.nv(j);i++)
      {
        // add the parent data for this vertex
        //parent.push_back( k ); // TODO
      }
    }
  }
}

template<typename type>
void
Topology<type>::get_triangles( Topology<Simplex>& triangulation ) const
{
  avro_assert( triangulation.points().nb() == 0 );
  avro_assert( triangulation.number() == 2 );
  avro_assert( triangulation.nb() == 0 );

  // the reason we add all the vertices is in case the caller doesn't want to
  // actually change the points
  for (index_t k=0;k<points_.nb();k++)
    triangulation.points().create( points_[k] );

  // setting the number to 2 will instruct the triangulator to only save
  // simplices once we get to triangles
  triangulate( triangulation );
}
#endif

template<typename type>
void
Topology<type>::facet( const index_t k , const index_t j ,
                       std::vector<index_t>& f ) const
{
  element_.facet( operator()(k) , j , f );
}

template<typename type>
void
Topology<type>::get_boundary( Topology<type>& boundary ) const
{
  avro_implement;
}

template<typename type>
void
Topology<type>::get_elements( Topology<type>& topology ) const
{
  std::vector<const Topology<type>*> children;
  this->get_children(children);
  children.push_back(this);
  for (index_t i=0;i<children.size();i++)
  {
    if (children[i]->number()!=topology.number()) continue;
    for (index_t k=0;k<children[i]->nb();k++)
      topology.add( (*children[i])(k) , children[i]->nv(k) );
  }
}

template<typename type>
bool
Topology<type>::ghost( const index_t* v , index_t nv ) const
{
  // TODO: make sure that ghost elements always have the first entry
  // as the ghost vertex even when the mesh is oriented
  for (index_t j=0;j<nv;j++)
  {
    if (v[j]<points_.nb_ghost())
      return true;
  }
  return false;
}

template<typename type>
bool
Topology<type>::ghost( index_t k ) const
{
  avro_assert_msg( k < nb() , "requested cell %lu but nb = %lu", k , nb() );
  return ghost( (*this)(k) , nv(k) );
}

template<typename type>
index_t
Topology<type>::nb_ghost() const
{
  index_t count = 0;
  for (index_t k=0;k<nb();k++)
  {
    if (ghost(k)) count++;
  }
  return count;
}

template<typename type>
index_t
Topology<type>::nb_real() const
{
  index_t count = 0;
  for (index_t k=0;k<nb();k++)
  {
    if (ghost(k)) continue;
    count++;
  }
  return count;
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
bool
Topology<type>::has( const index_t value ) const
{
  for (index_t k=0;k<data_.size();k++)
    if (data_[k]==value) return true;
  return false;
}

template<typename type>
index_t
Topology<type>::cardinality( const index_t* v , index_t nv ) const
{
  std::vector<index_t> d(v,v+nv);
  std::sort(d.begin(),d.end());
  index_t count = 0;
  bool has;
  std::vector<index_t> tk(nv);
  for (index_t k=0;k<nb();k++)
  {
    has = true;
    for (index_t j=0;j<nv;j++)
      tk[j] = (*this)(k,j);
    if (nv!=tk.size()) return false;
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
    avro_assert( number_+1 == points_.dim() );

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
    if (ghost(k)) continue;
    v += element_.volume( points_ , operator()(k) , nv(k) );
  }
  return v;
}

template<typename type>
void
Topology<type>::get_volumes( std::vector<real_t>& volumes ) const
{
  if (volumes.size()!=nb())
    volumes.resize( nb() );

  for (index_t k=0;k<nb();k++)
  {
    if (ghost(k))
      volumes[k] = 0.0;
    else
      volumes[k] = element_.volume( points_ , (*this)(k) , nv(k) );
  }
}

template<typename type>
void
Topology<type>::remove_point( const index_t k )
{
  // is this the last point?
  //bool last = ( k >= (points_.nb()-1) );
  points_.remove(k);

  //if (last) avro_implement;

  // decrement any indices higher than the original index
  for (index_t i=0;i<data_.size();i++)
  {
    if (data_[i]>k)
      data_[i]--;
  }
}

template<typename type>
void
Topology<type>::remove_elements( const std::vector<index_t>& elems0 )
{
  // make a copy so we can sort
  std::vector<index_t> elems(elems0.begin(),elems0.end());
  std::sort( elems.begin() , elems.end() );

  #if 0
  for (index_t k=0;k<elems.size();k++)
    remove( elems[k]-k );
  #elif 1
  std::vector<index_t> data0;
  //index_t n = 0;
  std::set<index_t> selem( elems.begin() , elems.end() );
  data0.reserve( data_.size() );
  for (index_t k = 0; k < nb(); k++) {
    if (selem.find(k) != selem.end() )
      continue;
    for (index_t j = 0 ; j < nv(k); j++)
      data0.push_back( (*this)(k,j) );
  }
  data0.shrink_to_fit();
  //avro_assert_msg( n == (number()+1)*(nb() - elems.size()) , "n = %lu, ndata = %lu" , n , (number()+1)*(nb()-elems.size()) );
  data_.assign( data0.begin() , data0.end() );
  #else
  std::list<index_t> ldata( data_.begin() , data_.end() );
  std::list<index_t>::iterator it = ldata.begin();
  for (index_t k = 0; k < elems.size(); k++) {
    advance(it,(number()+1)*(elems[k] -k));
    for (index_t j = 0; j < number()+1; j++) {
      advance(it,1);
      ldata.erase(it);
    }
    it = ldata.begin();
  }
  data_.assign( ldata.begin() , ldata.end() );

  #endif
}

template<typename type>
void
Topology<type>::remove_points( const std::vector<index_t>& pts0 )
{
  std::vector<index_t> pts(pts0.begin(),pts0.end());
  std::sort(pts.begin(),pts.end());

  //printf("before\n");
  //Table<index_t>::print();

  for (index_t i=0;i<pts.size();i++)
  {
    remove_point(0);
    //points_.remove( 0 );
  }

  for (index_t j=0;j<data_.size();j++)
  {
    //data_[j] -= pts.size();
  }

  //printf("after\n");
  //Table<index_t>::print();
  //points_.print();
}

template<typename type>
bool
Topology<type>::all_points_accounted() const
{
  std::vector<bool> accounted( points_.nb() );
  for (index_t k=0;k<data_.size();k++)\
    accounted[ data_[k] ] = true;

  bool result = true;
  for (index_t k=0;k<accounted.size();k++)
  {
    if (!accounted[k])
    {
      printf("point %lu is not accounted for\n",k);
      points_.print(k);
      result = false;
    }
  }
  return result;
}

template<typename type>
void
Topology<type>::determine_unused( std::vector<index_t>& pts ) const
{
  std::vector<bool> accounted( points_.nb() );
  for (index_t k=0;k<data_.size();k++)
    accounted[ data_[k] ] = true;

  for (index_t k=0;k<accounted.size();k++)
  {
    if (!accounted[k])
    {
      pts.push_back(k);
    }
  }
}

template<typename type>
void
Topology<type>::remove_unused( std::vector<index_t>* pidx0 )
{
  std::vector<bool> accounted( points_.nb() );
  for (index_t k=0;k<data_.size();k++)
    accounted[ data_[k] ] = true;

  std::vector<index_t> idx;
  std::vector<index_t>* pts;
  if (pidx0==nullptr)
    pts = &idx;
  else
    pts = pidx0;

  for (index_t k=0;k<accounted.size();k++)
  {
    if (!accounted[k])
    {
      pts->push_back(k);
    }
  }

#if 0
  for (index_t k=0;k<pts->size();k++)
    remove_point( (*pts)[k]-k);
#elif 1

  // this is inefficient
  move_to_front( *pts );

  // batch erase the points
  index_t n = pts->size();
  for (index_t j = 0; j < data_.size(); j++) {
    avro_assert( data_[j] >= n );
    data_[j] -= n;
  }
  points_.batch_erase(n);
#else

#endif
}

template<typename type>
void
Topology<type>::move_to_front( const std::vector<index_t>& pts , std::map<index_t,index_t>* pidx0 )
{
  std::map<index_t,index_t> idx;
  std::map<index_t,index_t>* pidx;
  if (pidx0==nullptr)
    pidx = &idx;
  else
    pidx = pidx0;

  for (index_t k=0;k<pts.size();k++)
  {
    pidx->insert( {pts[k],k} );
  }

  for (index_t j=0;j<data_.size();j++)
  {
    index_t offset = 0;
    for (index_t i=0;i<pts.size();i++)
    {
      if (data_[j]<pts[i]) offset++;
    }
    if (pidx->find(data_[j])==pidx->end())
      pidx->insert( {data_[j] , data_[j]+offset} );
  }

  // remap all the indices
  for (index_t j=0;j<data_.size();j++)
    data_[j] = pidx->at(data_[j]);

  // now move all the points to the front
  #if 0
  for (index_t k=0;k<pts.size();k++)
  {
    points_.move_to( pts[k] , k );
  }
  #else
  points_.move_to_front( pts );
  #endif

}

template<typename type>
void
Topology<type>::close()
{
  if (closed_)
  {
    printf("mesh is already closed!\n");
    return;
  }

  // compute the boundary
  Topology<type> bnd(this->points_,this->number_-1);
  get_boundary(bnd);

  // count how many bodies there are
  std::vector<int> bodies;
  const std::vector<index_t>& idx = bnd.data();
  for (index_t k=0;k<idx.size();k++)
  {
    if (points_.body( idx[k] )==0) continue; // interior
    if (points_.body( idx[k] )<0) continue; // partition boundary
    bodies.push_back( this->points_.body( idx[k] ) );
  }
  uniquify(bodies);

  // special case when nothing is tagged
  if (bodies.size()==0)
  {
    // tag all boundary vertices with some body (1)
    bodies.push_back(1);
    for (index_t k=0;k<bnd.nb();k++)
    for (index_t i=0;i<bnd.nv(k);i++)
      points_.body( bnd(k,i) ) = 1;
  }

  index_t nb_offset = 0;
  for (index_t i=0;i<bodies.size();i++)
  {
    if (bodies[i]<0) continue;
    nb_offset++;
    index_t id = this->points_.nb_ghost();
    this->points_.create_ghost();
    this->points_.body( id ) = bodies[i];
  }

  bnd.offset_by(nb_offset);
  offset_by(nb_offset);

  for (index_t k=0;k<nb();k++)
  for (index_t j=0;j<nv(k);j++)
    avro_assert( (*this)(k,j) < points_.nb() );

  // pick the ghost for partition boundaries
  int pbody = -1;
  for (index_t k=0;k<bodies.size();k++)
  {
    if (bodies[k]<0) continue;
    pbody = bodies[k];
    break;
  }
  avro_assert( pbody > 0 );

  // create ghost elements
  for (index_t k=0;k<bnd.nb();k++)
  {
    // look ahead which boundary this is
    avro_assert( bnd.nv(k)>0 );
    int ibnd = this->points().body( bnd(k,0) );

    // reset the body for partition boundaries
    if (ibnd<0) ibnd = pbody;
    ibnd = 1; // philip: july 15 -> bnd tag not even used

    std::vector<index_t> fk = bnd.get(k);
    std::sort(fk.begin(),fk.end());
    fk.insert( fk.begin() , index_t(ibnd-1) );
    this->add( fk.data() , fk.size() );
  }

  for (index_t k=0;k<nb();k++)
  for (index_t j=0;j<nv(k);j++)
    avro_assert( (*this)(k,j) < points_.nb() );

  closed_ = true;
}

template<typename type>
void
Topology<type>::intersect( const std::vector<index_t>& facet , std::vector<index_t>& elems ) const
{
  if (facet.size()==1)
    inverse_.ball( facet[0] , elems );
  else if (facet.size()==2)
    inverse_.shell( facet[0] , facet[1] , elems );
  else if (facet.size()==3)
    inverse_.shell( facet[0] , facet[1] , facet[2] , elems );
  else
    avro_implement;
}

template<typename type>
void
Topology<type>::set_points( Points& points )
{
  points_ = points;
  //points_.print();
  for (index_t k=0;k<nb_children();k++)
    Tree<Topology<type>>::child(k).set_points(points);
}

template<typename type>
void
Topology<type>::print_header() const
{
  printf("topology %p: order = %u, number = %u\n",(void*)this,element_.order(),element_.number());
}

template<typename type>
void
Topology<type>::build_structures()
{
  neighbours_.fromscratch() = true;
  neighbours_.forceCompute();
  neighbours_.compute();
  inverse_.build();
}

} // avro
