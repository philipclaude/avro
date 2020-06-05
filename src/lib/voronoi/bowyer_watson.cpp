//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "adaptation/cavity.h"

#include "library/ckf.h"

#include "mesh/points.h"
#include "mesh/search.h"

#include "voronoi/algorithms.h"

#include <triangle/predicates.h>

namespace avro
{

BowyerWatson::BowyerWatson( Points& delaunay ) :
  Topology<Simplex>( points_ , delaunay.dim() ),
  points_( delaunay.dim() ),
  delaunay_(delaunay),
  cloud_(delaunay)
{}

void
BowyerWatson::initialize()
{
  kdtree_ = initializeKdTree(cloud_);
  kdtree_->build();

  std::vector<real_t> xmin( points_.dim() ,  1e20 );
  std::vector<real_t> xmax( points_.dim() , -1e20);

  for (index_t k=0;k<delaunay_.nb();k++)
  {
    for (coord_t d=0;d<delaunay_.dim();d++)
    {
      real_t xd = delaunay_[k][d];
      if (xd < xmin[d]) xmin[d] = xd;
      if (xd > xmax[d]) xmax[d] = xd;
    }
  }

  std::vector<real_t> xmid( points_.dim() );
  for (coord_t d=0;d<points_.dim();d++)
  {
    xmin[d] = xmin[d]*1.01 -0.01;
    xmax[d] = xmax[d]*1.01 +0.01;
    xmid[d] = xmin[d] + 0.5*( xmax[d] - xmin[d] );
  }

  std::vector<index_t> dims( points_.dim() , 2 );
  CKF_Triangulation cube(dims);

  // center the cube at the midpoint (they are initially from 0 -> 1)
  for (index_t k=0;k<cube.points().nb();k++)
  {
    for (index_t d=0;d<cube.points().dim();d++)
    {
      cube.points()[k][d] = xmin[d] + (xmax[d] - xmin[d])*cube.points()[k][d];
    }
    fake_.push_back( points_.nb() );
    points_.create( cube.points()[k] );
  }

  for (index_t k=0;k<cube.nb();k++)
    add( cube(k) , cube.nv(k) );

  close();
  orient();
  this->neighbours().compute();
  this->inverse().build();
}

real_t
BowyerWatson::insphere( index_t elem , index_t point )
{
  const index_t* t = operator()(elem);
  if (points_.dim()==2)
  {
    return ::incircle( points_[t[0]] , points_[t[1]] , points_[t[2]] , points_[point] );
  }
  else if (points_.dim()==3)
  {
    return -::insphere( points_[t[0]] , points_[t[1]] , points_[t[2]] , points_[t[3]], points_[point] );
  }
  else if (points_.dim()==4)
  {
    avro_implement;
    //return -::inhypersphere( points_[t[0]] , points_[t[1]] , points_[t[2]] , points_[t[3]], , points_[t[4]] , points_[point] );
  }
  else
    avro_implement;
}

void
BowyerWatson::search( index_t point , index_t elem , std::set<index_t>& elems )
{
  for (index_t j=0;j<index_t(number_+1);j++)
  {
    int neighbour = neighbours()(elem,j);
    if (neighbour<0) continue;
    if (ghost(neighbour)) continue;
    if (elems.find(neighbour)!=elems.end()) continue;

    if (insphere(neighbour,point)>0)
    {
      elems.insert( neighbour );
      search(point,neighbour,elems);
    }
  }
}

void
BowyerWatson::compute()
{
  // compute super-triangle
  initialize();

  // create a cavity operator to handle the removal/insertion
  Cavity<Simplex> cavity( *this );
  cavity.set_enlarge(false);
  cavity.check_visibility(false);

  // build a kdtree so we insert vertices close to the last one
  std::vector<index_t> idx( delaunay_.nb() );
  std::vector<real_t> dist( delaunay_.nb() );
  index_t nu = delaunay_.nb();
  kdtree_->getNearestNeighbours( delaunay_[0] , idx , dist , nu );

  ElementSearch<Simplex> searcher(*this);

  // add every points one by one
  for (index_t k=0;k<delaunay_.nb();k++)
  {

    cavity.clear();

    // add the point into the triangulation
    index_t point = points_.nb();
    points_.create( delaunay_[idx[k]] );
    inverse().create(1);

    // find all triangles whose circumsphere contains this point
    std::set<index_t> elems;
    index_t start = nb();
    if (k==0) // first point
    {
      for (index_t j=0;j<nb();j++)
      {
        if (ghost(j)) continue;
        if (insphere(j, point )>0)
        {
          start = j;
          break;
        }
      }
    }
    else
    {
      // get the ball of the previous vertex since this is close to the one
      // currently being added
      std::vector<index_t> ball;
      inverse().ball( point-1 , ball );

      // check if any simplex in the ball of the previous vertex encloses
      // the current point
      for (index_t j=0;j<ball.size();j++)
      {
        if (ghost(ball[j])) continue;
        if (insphere(ball[j], point )>0)
        {
          start = ball[j];
          break;
        }
      }
      if (start==nb())
      {
        // the ball of the previous point does not enclose the current point
        // so search the mesh for an element enclosing the point
        int result = searcher.find( points_[point] , ball[0] );
        avro_assert( result>=0 );
        start = result;
      }
    }
    avro_assert( start < nb() );

    // search for all elements broken by this point, starting with the one we found
    elems.insert(start);
    search( point , start , elems );

    // apply the cavity operator to this point
    // (connect the point to the boundary and insert into triangulation)
    std::vector<index_t> C(elems.begin(),elems.end());
    bool result = cavity.compute( point , points_[point] , C );
    avro_assert( result );

    // apply the change to the topology
    apply(cavity);
  }

  // delete the points added from the initial supertriangle
  // as well as any triangle connected to one of these fake points
  std::vector<index_t> removed;
  for (index_t k=0;k<nb();k++)
  {
    for (index_t j=0;j<fake_.size();j++)
      if (has(k,fake_[j]))
        removed.push_back(k);
  }
  uniquify(removed);

  for (index_t k=0;k<removed.size();k++)
    remove( removed[k] - k );

  for (index_t k=0;k<fake_.size();k++)
    remove_point( fake_[k] - k );
}

} // avro
