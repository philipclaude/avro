//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_GRAPHICS_PRIMITIVE_H_
#define avro_LIB_GRAPHICS_PRIMITIVE_H_

#include "common/tree.h"

#include "graphics/gl.h"

#include <map>

namespace avro
{

template<typename Master_t> class Topology;
class TopologyBase;

namespace graphics
{

class ShaderProgram;
class Plotter;
class GraphicsManager;
class SceneGraph;
class ClippingPlane;

class Primitive : public Tree<Primitive>
{
public:
  Primitive( const TopologyBase& topology , SceneGraph* scene );

  coord_t number() const { return number_; }

  void write( GraphicsManager& manager , const ClippingPlane* plane );
  void extract(const ClippingPlane* plane=nullptr);

  bool& visible() { return visible_; }
  bool& points_on() { return points_on_; }
  bool& edges_on() { return edges_on_; }
  bool& triangles_on() { return triangles_on_; }
  float& transparency() { return transparency_; }

  const std::vector<index_t>& triangles() const { return triangles_; }
  const std::vector<index_t>& edges() const { return edges_; }
  const std::vector<real_t>& points() const { return points_; }
  const std::vector<real_t>& colors() const { return colors_; }
  const std::vector<real_t>& normals() const { return normals_; }

  SceneGraph* scene() { return scene_; }
  const TopologyBase& topology() const;

  void hide(bool hidden);
  bool& hidden() { return hidden_; }
  void hide();
  void show();
  void set_active( const std::string& x , index_t rank=0 );

  void get_field_limits( real_t* lims ) const;

protected:
  coord_t number_;
  const TopologyBase& topology_;

  index_t rank_;
  std::string active_;
  SceneGraph* scene_;
  bool visible_;

  std::vector<index_t> edges_;
  std::vector<index_t> triangles_;
  std::vector<real_t>  points_;
  std::vector<real_t>  normals_;
  std::vector<real_t>  colors_;

  bool triangles_on_;
  bool edges_on_;
  bool points_on_;
  float transparency_;

  bool hidden_;

  real_t ulim_[2];
};

} // graphics

} // avro

#endif
