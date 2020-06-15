//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "common/tree.hpp"

#include "graphics/application.h"
#include "graphics/colormap.h"
#include "graphics/primitive.h"

#include "mesh/decomposition.h"
#include "mesh/field.h"
#include "mesh/points.h"
#include "mesh/topology.h"

namespace avro
{

namespace graphics
{

Primitive::Primitive( const TopologyBase& topology , SceneGraph* scene ) :
  number_(topology.number()),
  topology_(topology),
  rank_(0),
  active_(""),
  scene_(scene),
  visible_(true),
  triangles_on_(true),
  edges_on_(true),
  points_on_(false),
  transparency_(1.0),
  hidden_(false)
{}

const TopologyBase&
Primitive::topology() const
{
  return topology_;
}

void
Primitive::extract( const ClippingPlane* plane )
{
  index_t nb_triangles = 0;
  coord_t dim0 = topology_.points().dim();
  coord_t number = topology_.number();

  points_.clear();
  edges_.clear();
  triangles_.clear();
  normals_.clear();
  colors_.clear();

  const real_t* center = scene_->focus();
  const real_t scale = 1./center[3];

  if (number==0)
  {
    // get the edges from the topology
    for (index_t k=0;k<topology_.nb();k++)
    {
      for (index_t j=0;j<dim0;j++)
        points_.push_back( scale*( topology_.points()[topology_(k,0)][j] -center[j]) );
      for (index_t j=dim0;j<3;j++)
        points_.push_back( scale*( 0.0 - center[j] ) );

      points_on_ = true;
    }
    return;
  }

  if (number==1)
  {
    // get the edges from the topology
    std::vector<index_t> edges;
    topology_.get_edges( edges , plane );

    std::vector<index_t> point_map0( topology_.points().nb() );
    for (index_t k=0;k<edges.size();k++)
    {
      for (index_t j=0;j<dim0;j++)
        points_.push_back( scale*( topology_.points()[edges[k]][j]  - center[j]) );
      for (index_t j=dim0;j<3;j++)
        points_.push_back( scale*( 0.0 - center[j]) );

      edges_.push_back( edges_.size() );
    }
    return;
  }

  // get the triangles from the topology
  std::shared_ptr<SimplicialDecompositionBase> pdecomposition;
  if (topology_.type_name()=="simplex")
    pdecomposition = std::make_shared<SimplicialDecomposition<Simplex>>(*dynamic_cast<const Topology<Simplex>*>(&topology_));
  else
    pdecomposition = std::make_shared<SimplicialDecomposition<Polytope>>(*dynamic_cast<const Topology<Polytope>*>(&topology_));
  pdecomposition->extract(plane);
  const SimplicialDecompositionBase& decomposition = *pdecomposition.get();

  const Points& points = decomposition.points();

  // TODO only account for triangles that are visible to the clipping plane
  std::vector<index_t> triangles;
  std::vector<index_t> parents;
  decomposition.get_simplices(2,triangles,parents); // setting 2 retrieves triangles
  index_t nb_points = triangles.size(); // duplicated points for cell-based colors

  // allocate the vertex data
  colors_.resize( 3*nb_points, 0 );
  normals_.resize( 3*nb_points, 0 );

  std::vector<index_t> point_map0( points.nb() );
  std::vector<index_t> point_map1;
  for (index_t k=0;k<triangles.size();k++)
  {
    for (index_t j=0;j<dim0;j++)
      points_.push_back( scale*(points[triangles[k]][j] - center[j]) );
    for (index_t j=dim0;j<3;j++)
      points_.push_back( scale*(0.0                     - center[j]) );

    point_map0[ triangles[k] ] = triangles_.size();
    point_map1.push_back( triangles[k] );

    if (number>=2)
      triangles_.push_back( triangles_.size() );
  }
  nb_triangles = triangles_.size()/3;

  if (number>=1)
  {
    // get the edges from the topology
    std::vector<index_t> edges;
    topology_.get_edges( edges , plane );
    for (index_t k=0;k<edges.size();k++)
      edges_.push_back( point_map0[ edges[k] ] );
  }

  std::vector<real_t> U;
  real_t umin,umax;

  bool constant_color = false;
  float color[3];
  const Fields& fields = topology_.fields();
  if (fields.has(active_))
  {
    // retrieve the active field
    umin = fields[active_].min(rank_);
    umax = fields[active_].max(rank_);

    // evaluate the active field (with rank) on the triangulation points
    const std::vector<index_t>& point_parents = decomposition.point_parents();
    const Table<real_t>& reference_coordinates = decomposition.reference_coordinates();
    fields[active_].evaluate( rank_ , point_parents , reference_coordinates , U );
  }
  else
  {
    // constant color
    umin = 0.0;
    umax = 1.0;
    U.resize( decomposition.points().nb() , 1.0 );
    color[0] = 0.5;
    color[1] = 0.5;
    color[2] = 0.5;
    constant_color = true;
  }
  float lims[2] = {float(umin),float(umax)};
  ulim_[0] = umin;
  ulim_[1] = umax;

  // compute the color of the (unduplicated) triangulation points
  Colormap& colormap = scene_->colormap();
  colormap.set_limits(lims);
  std::vector<index_t> color0( points.nb()*3 );
  for (index_t k=0;k<triangles.size()/3;k++)
  {
    // retrieve the parent
    for (index_t j=0;j<3;j++)
    {
      index_t p = triangles[3*k+j];
      real_t  u = U[p];
      if (!constant_color)
        colormap.map( u , color );

      color0[3*p  ] = 255*color[0];
      color0[3*p+1] = 255*color[1];
      color0[3*p+2] = 255*color[2];
    }
  }

  // map the colors to all the vertices
  for (index_t k=0;k<nb_points;k++)
  {
    for (index_t j=0;j<3;j++)
      colors_[3*k+j] = color0[3*point_map1[k]+j];
  }

  // compute the normals
  if (topology_.number()>=2)
  {
    real_t x1,y1,z1,x2,y2,z2,x3,y3,z3;
    real_t dis;
    real_t xnor,ynor,znor;
    std::vector<index_t> count( nb_points, 0 );
    for (index_t k=0;k<nb_triangles;k++)
    {
      x1 = points_[ 3*triangles_[3*k  ] + 0 ];
      y1 = points_[ 3*triangles_[3*k  ] + 1 ];
      z1 = points_[ 3*triangles_[3*k  ] + 2 ];
      x2 = points_[ 3*triangles_[3*k+1] + 0 ];
      y2 = points_[ 3*triangles_[3*k+1] + 1 ];
      z2 = points_[ 3*triangles_[3*k+1] + 2 ];
      x3 = points_[ 3*triangles_[3*k+2] + 0 ];
      y3 = points_[ 3*triangles_[3*k+2] + 1 ];
      z3 = points_[ 3*triangles_[3*k+2] + 2 ];

      xnor = (y3-y1)*(z2-z1)-(y2-y1)*(z3-z1);
      ynor = (z3-z1)*(x2-x1)-(x3-x1)*(z2-z1);
      znor = (x3-x1)*(y2-y1)-(y3-y1)*(x2-x1);

      dis  = sqrtf(xnor*xnor + ynor*ynor + znor*znor);

      if (dis != 0.0f)
      {
        normals_[ 3*triangles_[3*k  ] + 0 ] += xnor/dis;
        normals_[ 3*triangles_[3*k  ] + 1 ] += ynor/dis;
        normals_[ 3*triangles_[3*k  ] + 2 ] += znor/dis;
        normals_[ 3*triangles_[3*k+1] + 0 ] += xnor/dis;
        normals_[ 3*triangles_[3*k+1] + 1 ] += ynor/dis;
        normals_[ 3*triangles_[3*k+1] + 2 ] += znor/dis;
        normals_[ 3*triangles_[3*k+2] + 0 ] += xnor/dis;
        normals_[ 3*triangles_[3*k+2] + 1 ] += ynor/dis;
        normals_[ 3*triangles_[3*k+2] + 2 ] += znor/dis;

        count[ triangles_[3*k  ] ]++;
        count[ triangles_[3*k+1] ]++;
        count[ triangles_[3*k+2] ]++;
      }
    }

    // normalize again
    for (index_t k=0;k<nb_points;k++)
    {
      if (count[k] <= 1) continue;
      dis  = count[k];
      xnor = normals_[3*k  ] / dis;
      ynor = normals_[3*k+1] / dis;
      znor = normals_[3*k+2] / dis;
      dis  = sqrtf(xnor*xnor + ynor*ynor + znor*znor);
      normals_[3*k  ] = xnor/dis;
      normals_[3*k+1] = ynor/dis;
      normals_[3*k+2] = znor/dis;
    }
  }
  else if (topology_.number()==1)
  {
    printf("setting color for lines!!\n");
    colors_.resize( 3*nb_points, 255.0f );
  }
}

void
Primitive::hide(bool hidden)
{
  triangles_on_ = !hidden;
  edges_on_ = !hidden;

  for (index_t k=0;k<nb_children();k++)
    child(k).hide(hidden);
}

void
Primitive::hide()
{
  hidden_ = true;
  hide(true);
}

void
Primitive::show()
{
  hidden_ = false;
  hide(false);
}

void
Primitive::write( GraphicsManager& manager , const ClippingPlane* plane )
{
  extract(plane);
  manager.write( *this );

  for (index_t k=0;k<nb_children();k++)
    child(k).write(manager,plane);
}

void
Primitive::set_active( const std::string& x , index_t rank )
{
  active_ = x;
  rank_   = rank;
  for (index_t k=0;k<nb_children();k++)
    child(k).set_active(x,rank);
}

void
Primitive::get_field_limits( real_t* lims ) const
{
  lims[0] = std::min( lims[0] , ulim_[0] );
  lims[1] = std::max( lims[1] , ulim_[1] );
  for (index_t k=0;k<nb_children();k++)
    child(k).get_field_limits( lims );
}


} // graphics

template class Tree<graphics::Primitive>;

} // avro
