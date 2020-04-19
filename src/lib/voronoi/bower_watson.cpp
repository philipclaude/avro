#include "adaptation/cavity.h"

#include "library/ckf.h"

#include "mesh/points.h"

#include "voronoi/algorithms.h"

#include <triangle/predicates.h>

namespace avro
{

BowerWatson::BowerWatson( Points& delaunay ) :
  Topology<Simplex>( points_ , delaunay.dim() ),
  points_( delaunay.dim() ),
  delaunay_(delaunay)
{}

void
BowerWatson::initialize()
{
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
BowerWatson::insphere( index_t elem , index_t point )
{
  const index_t* t = operator()(elem);
  if (points_.dim()==2)
  {
    return ::incircle( points_[t[0]] , points_[t[1]] , points_[t[2]] , points_[point] );
  }
  if (points_.dim()==3)
  {
    return -::insphere( points_[t[0]] , points_[t[1]] , points_[t[2]] , points_[t[3]], points_[point] );
  }
  else
    avro_implement;
}

void
BowerWatson::compute()
{
  // compute super-triangle
  initialize();

  // create a cavity operator to handle the removal/insertion
  Cavity<Simplex> cavity( *this );
  cavity.set_enlarge(false);
  cavity.check_visibility(false);

  // add every points one by one
  for (index_t k=0;k<delaunay_.nb();k++)
  {
    cavity.clear();

    // add the point into the triangulation
    index_t point = points_.nb();
    points_.create( delaunay_[k] );
    inverse().create(1);

    // find all triangles whose circumsphere contains this point
    std::vector<index_t> elems;

    // TODO more efficient search through mesh
    // using the simplex connectivity (since this is currently O(n^2))
    for (index_t j=0;j<nb();j++)
    {
      if (ghost(j)) continue;
      if (insphere(j, point )>0)
        elems.push_back(j);
    }
    avro_assert( elems.size()>0 );

    // apply the cavity operator to this point
    // (connect the point to the boundary and insert into triangulation)
    bool result = cavity.compute( point , points_[point] , elems );
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
