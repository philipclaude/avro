//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "common/kdtree.h"

#include "mesh/points.h"

#include <nanoflann/nanoflann.hpp>

#include <vector>

namespace avro
{

PointCloud::PointCloud( Points& points ) :
  points_(points)
{}

real_t
PointCloud::kdtree_get_pt(const size_t idx, int d) const
 { return points_[idx][d]; }

size_t
PointCloud::kdtree_get_point_count() const
 { return points_.nb(); }

coord_t
PointCloud::dim() const
  { return points_.dim(); }

template<coord_t dim>
KdTree<dim>::KdTree( PointCloud& _cloud ) :
  KdTreeNd(_cloud)
{}

template<coord_t dim>
void
KdTree<dim>::build()
{
  tree_ = std::make_shared<nanoflann::KDTreeSingleIndexAdaptor<
                nanoflann::L2_Simple_Adaptor<real_t,PointCloud>,
                PointCloud,dim> >
                (dim,cloud_,nanoflann::KDTreeSingleIndexAdaptorParams(10));
  tree_->buildIndex();
}

template<coord_t dim>
void
KdTree<dim>::getNearestNeighbours( const real_t* q, std::vector<index_t>& idx ,
                                   std::vector<real_t>& dist2 , index_t& nu )
{
  avro_assert( idx.size()==dist2.size() );
  nu = tree_->knnSearch(q,idx.size(),idx.data(),dist2.data());
  //idx.resize(nu);
  //dist2.resize(nu);
}


template class KdTree<1>;
template class KdTree<2>;
template class KdTree<3>;
template class KdTree<4>;
template class KdTree<5>;
template class KdTree<6>;
template class KdTree<7>;
template class KdTree<8>;
template class KdTree<9>;
template class KdTree<10>;
template class KdTree<11>;
template class KdTree<12>;
template class KdTree<13>;
template class KdTree<14>;

} // avro
