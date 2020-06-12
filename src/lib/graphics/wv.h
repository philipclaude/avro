//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_GRAPHICS_WV_H_
#define avro_LIB_GRAPHICS_WV_H_

#include "common/types.h"

#include "graphics/manager.h"

#include <wsss.h>

#include <map>

#ifdef __cplusplus
extern "C"
{
#endif
void browserMessage( void* wsi , char* text , /*@unused@*/ int lena );
#ifdef __cplusplus
}
#endif

namespace avro
{

namespace graphics
{

class WV_Manager : public GraphicsManager
{
  template<typename type> friend class Application;

private:
  WV_Manager();

  void write( Primitive& primitive );
  void write( const std::string& name , coord_t number , const std::vector<real_t>& points , const std::vector<index_t>& edges , const std::vector<index_t>& triangles , const std::vector<real_t>& colors )
  { avro_assert_not_reached; }
  void draw( SceneGraph& scene , TransformFeedbackResult* feedback=nullptr )
  { avro_assert_not_reached; }
  void draw( const std::string& name , coord_t number , const DrawingParameters& params )
  { avro_assert_not_reached; }

  void select_shader( const std::string& name , const std::string& shader_name )
  { avro_assert_not_reached; }

  wvContext* context() { return context_; }

private:

  wvContext* context_;

  // map from primitive to wv gprim index
  std::map<Primitive*,index_t> index_;

  int bias_;
  float fov_,znear_,zfar_;
  float eye_[3],center_[3],up_[3];


};

} // graphics

} // avro

#endif
