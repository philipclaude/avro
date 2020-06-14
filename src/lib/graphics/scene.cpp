//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "graphics/controls.h"
#include "graphics/scene.h"

#include "mesh/points.h"

#include <limits>

namespace avro
{

namespace graphics
{

SceneGraph::SceneGraph() :
  update_(true)
{
  menu_["primitives"] = {};
  focus_[0] = focus_[1] = focus_[2] = 0.0;
  focus_[3] = 1.0;
}

index_t
SceneGraph::add_primitive( const TopologyBase& topology )
{
  // create a primitive for the incoming topology
  index_t id = primitive_.size();
  Primitive_ptr primitive = std::make_shared<Primitive>(topology,this);
  primitive_.push_back(primitive);

  // create graphics primitives for the topology children too!
  std::vector<const TopologyBase*> children;
  topology.get_topologies(children);
  for (index_t k=0;k<children.size();k++)
  {
    primitive->add_child( std::make_shared<Primitive>(*children[k],this) );
  }
  return id;
}

void
SceneGraph::remove( index_t k )
{
  // warning: i don't think this is tested
  primitive_.erase( primitive_.begin()+k );
}

void
SceneGraph::write( GraphicsManager& manager , const ClippingPlane* plane  )
{
  for (index_t k=0;k<primitive_.size();k++)
    primitive_[k]->write(manager,plane);

  std::vector<std::string> primitives;
  for (index_t k=0;k<primitive_.size();k++)
  {
    MenuEntry entry(*primitive_[k].get());
    primitives.push_back( entry.dump() );
  }
  menu_["primitives"] = primitives;
}

void
SceneGraph::update_matrices( const Controls& controls )
{
  normal_matrix_ = controls.normal();
  mvp_matrix_ = controls.model_view_projection();
}

void
SceneGraph::get_bounding_box( real_t* box ) const
{
  for (index_t k=0;k<primitive_.size();k++)
  {
    const Points& points = primitive_[k]->topology().points();
    coord_t dim = ( points.dim() < 3 ) ? points.dim() : 3;
    for (index_t j=0;j<points.nb();j++)
    {
      for (coord_t d=0;d<dim;d++)
      {
        box[  d] = std::min( box[d  ] , points[j][d] );
        box[3+d] = std::max( box[3+d] , points[j][d] );
      }
    }
    for (coord_t d=dim;d<3;d++)
    {
      box[ d]  = std::min( box[d  ] , -1.0 );
      box[3+d] = std::max( box[3+d] , +1.0 );
    }
  }
}

void
SceneGraph::get_color_limits( real_t* clim ) const
{
  printf("nb prims  = %lu\n",primitive_.size());
  for (index_t k=0;k<primitive_.size();k++)
  {
    printf("processing primitive %lu\n",k);
    std::vector<const Primitive*> children;
    primitive_[k]->get_children(children);
    for (index_t j=0;j<children.size();j++)
    {
      printf("getting children %lu\n",j);
      const std::vector<real_t>& color = primitive_[j]->colors();
      if (color.size()==0) continue;
      clim[0] = std::min( clim[0] , *std::min(color.begin(),color.end()) );
      clim[1] = std::max( clim[1] , *std::max(color.begin(),color.end()) );
    }
  }

}

void
SceneGraph::set_focus( real_t* focus )
{
  for (coord_t d=0;d<4;d++)
    focus_[d] = focus[d];
}

} // graphics

} // avro
