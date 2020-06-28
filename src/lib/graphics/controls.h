//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef AVRO_LIB_GRAPHICS_CONTROLS_H_
#define AVRO_LIB_GRAPHICS_CONTROLS_H_

#include "graphics/gl.h"
#include "graphics/math.h"

namespace avro
{

namespace graphics
{

class Controls
{

public:

  Controls( float fov , int width , int height , float znear=0.1f , float zfar=100.0f );

  void reset();
  bool update();
  void calculate_view();

  void mouse_down(int button, int action, int mods,int xpos,int ypos);
  void mouse_up();
  void mouse_move(int x, int y);
  void key_down(int key);
  void mouse_wheel(double xoffset ,double yoffset);

  const mat4& model_view_projection() const {return model_view_projection_; }
  const mat4& model_view() const { return model_view_; }
  const mat4& ui_matrix() const { return ui_matrix_; }
  const mat4& normal() const { return normal_; }
  const mat4& perspective() const { return perspective_; }

  void set_ui_matrix( const mat4& m ) { ui_matrix_ = m; }
  void set_mv_matrix( const mat4& m ) { model_view_ = m; }

  bool dragging;
  int modifier;

  void disable() { enabled_ = false; }
  void enable() { enabled_ = true; }

  const mat4& rotation() const { return rotation_; }
  const mat4& translation() const { return translation_; }

private:
  mat4 perspective_;
  mat4 model_view_;
  mat4 model_view_projection_;
  mat4 normal_;
  mat4 ui_matrix_;

  int width_;
  int height_;

  int offleft_;
  int offtop_;

  int startx_;
  int starty_;

  int cursorx_;
  int cursory_;

  float scale_;

  const vec3 eye_ = {0,0,7};
  const vec3 center_ = {0,0,0};
  const vec3 up_ = {0,1,0};

  bool enabled_;

  // rigid-body transformations (useful for defining clipping plane manipulations)
  mat4 rotation_;
  mat4 translation_;
};

} // graphics

} // avro

#endif
