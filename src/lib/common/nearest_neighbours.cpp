//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "common/kdtree.h"
#include "common/nearest_neighbours.h"
#include "common/tools.h"

#include "mesh/points.h"

#include "numerics/geometry.h"

#include <algorithm>

namespace avro
{

NearestNeighbours::NearestNeighbours( Points& _points , const index_t _knear ) :
  points_(_points) , knear_(_knear),
  neighbours_(TableLayout_Jagged)
{
  if (_knear==0) knear_ = points_.nb();
  if (knear_>=points_.nb()) knear_ = points_.nb();
  //neighbours_.set_rank(knear_);
  compute();
}

void
NearestNeighbours::compute()
{
  neighbours_.clear();

  index_t nv = points_.nb();

  #if 0
  coord_t nd = points_.dim();

  for (index_t k=0;k<nv;k++)
  {
    // compute all the distances
    std::vector<real_t> d( nv , 0. );
    std::vector<index_t> nn( nv );
    for (index_t j=0;j<nv;j++)
    {
      d[j]  = numerics::distance2( points_[k] , points_[j] , nd );
      nn[j] = j;
    }

    // sort by distance
    std::sort( nn.begin() , nn.end() , SortBy<real_t>(d) );

    // add the result
    neighbours_.add( nn.data() , knear_ );
  }

  #else

  PointCloud cloud( points_ );
  std::shared_ptr<KdTreeNd> kdtree = initializeKdTree( cloud );
  kdtree->build();

  for (index_t k=0;k<nv;k++)
  {
    std::vector<real_t> d( knear_ , 0. );
    std::vector<index_t> nn( knear_ );

    kdtree->getNearestNeighbours( points_[k] , nn , d , knear_ );
    neighbours_.add( nn.data() , nn.size() );
  }

  #endif
}

void
NearestNeighbours::print() const
{
  for (index_t k=0;k<points_.nb();k++)
  {
    printf("neighbours[%d] = (",int(k));
    for (index_t j=0;j<knear_;j++)
      printf(" %d",int(neighbours_(k,j)));
    printf(" )\n");
  }
}

} // avro
