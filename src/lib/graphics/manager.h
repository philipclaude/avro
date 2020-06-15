//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_GRAPHICS_MANAGER_H_
#define avro_LIB_GRAPHICS_MANAGER_H_

#include "common/error.h"
#include "common/types.h"

#include "graphics/listener.h"
#include "graphics/math.h"

namespace avro
{

namespace graphics
{

class Primitive;
class SceneGraph;
class TransformFeedbackResult;

typedef struct
{
  float transparency;
  mat4  mvp;
} DrawingParameters;

class GraphicsManager
{
public:
  virtual ~GraphicsManager() {}
  virtual void write( Primitive& primitive ) = 0;
  virtual void write( const std::string& name , coord_t number , const std::vector<real_t>& points , const std::vector<index_t>& edges , const std::vector<index_t>& triangles , const std::vector<real_t>& colors ) = 0;
  virtual void draw( SceneGraph& scene , TransformFeedbackResult* feedback=nullptr ) = 0;
  virtual void draw( const std::string& name , coord_t number , const DrawingParameters& params ) = 0;
  virtual void remove( const std::string& name ) = 0;

  Listener& listener() { return listener_; }

  virtual void select_shader( const std::string& name , const std::string& shader_name ) = 0;

private:
  Listener listener_;
};

class Vulkan_Manager : public GraphicsManager
{
  template<typename type> friend class Application;

private:
  Vulkan_Manager() {}

  void write( Primitive& primitive ) { avro_implement; }
  void write( const std::string& name , coord_t number , const std::vector<real_t>& points , const std::vector<index_t>& edges , const std::vector<index_t>& triangles , const std::vector<real_t>& colors )
  { avro_assert_not_reached; }
  void draw( SceneGraph& scene , TransformFeedbackResult* feedback=nullptr )
  { avro_implement; }

  void draw( const std::string& name , coord_t number , const DrawingParameters& params )
  { avro_implement; }

  void remove( const std::string& ) { avro_assert_not_reached; }

  void create_shaders()
  { avro_implement; }

  void select_shader( const std::string& name , const std::string& shader_name )
  { avro_implement; }

};

} // graphics

} // avro

#endif
