//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_GRAPHICS_WINDOW_H_
#define avro_LIB_GRAPHICS_WINDOW_H_

#include "graphics/clipping.h"
#include "graphics/controls.h"
#include "graphics/gl.h"
#include "graphics/math.h"

#include <memory>
#include <string>
#include <vector>

namespace avro
{

class Colormap;

namespace graphics
{

class SceneGraph;
class Interface;

#ifdef AVRO_WITH_GL

class GLFW_Window
{
public:
  GLFW_Window( GraphicsManager& manager , int width , int height , const std::string& title );

  void setup();

  // controls/interface functions
  void mouse_button_callback(int button,int action,int mods);
  void mouse_move_callback(double xpos, double ypos);
  void mouse_scroll_callback(double xpos, double ypos);
  void key_callback(int key, int scancode, int action, int mods);
  void create_interface();

  void make_current();
  void update_view();
  bool should_close();

  void save_eps( const std::string& filename );

  bool& changed() { return updated_; }

  void begin_draw();
  void draw_interface() const;
  void write_axes();
  void draw_axes();
  void draw_plane(const real_t* focus);
  void draw_colorbar(const Colormap& colormap,const real_t* ulim);
  void end_draw();

  index_t nb_scene() const { return scene_.size(); }
  SceneGraph& scene( index_t k ) { return *scene_[k].get(); }
  const SceneGraph& scene( index_t k ) const { return *scene_[k].get(); }

  index_t create_scene();

  Interface& interface() { return *interface_.get(); }
  const Interface& interface() const { return *interface_.get(); }

  Controls& controls() { return controls_; }
  Controls& clip_controls() { return clip_controls_; }
  ClippingPlane& clip_plane() { return clip_plane_; }

  GLFWwindow* window() { return window_; }
  const GLFWwindow* window() const { return window_; }

  int width() const { return width_; }
  int height() const { return height_; }

  index_t fps() const { return fps_; }
  void set_fps( index_t fps ) { fps_ = fps; }

  bool& modify_clipping_plane();
  bool& show_clipping_plane() { return show_clip_plane_; }
  void flip_clipping_normal();

  bool& show_axes() { return show_axes_; }
  bool& center_axes() { return center_axes_; }

  void clip();
  void reset_clip();

private:
  std::string title_;
  GLFWwindow* window_;
  GLFWimage images_[2];
  GraphicsManager& manager_;

  int width_;
  int height_;
  float fov_ = 45.0f;

  vec3 position_;
  float angles_[2];

  std::vector<std::shared_ptr<SceneGraph>> scene_;
  Controls controls_;
  Controls clip_controls_;

  bool modify_clip_plane_;
  bool show_clip_plane_;

  bool show_axes_;
  bool center_axes_;

  ClippingPlane clip_plane_;

  std::shared_ptr<Interface> interface_;

  bool updated_;

  index_t fps_;
};

#endif

} // graphics

} // avro

#endif
