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

#include "graphics/controls.h"
#include "graphics/gl.h"
#include "graphics/math.h"

#include <memory>
#include <string>
#include <vector>

namespace avro
{

namespace graphics
{

class SceneGraph;
class Interface;

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

  void begin_draw();
  void draw_interface() const;
  void end_draw();

  index_t nb_scene() const { return scene_.size(); }
  SceneGraph& scene( index_t k ) { return *scene_[k].get(); }
  const SceneGraph& scene( index_t k ) const { return *scene_[k].get(); }

  index_t create_scene();

  Interface& interface() { return *interface_.get(); }
  const Interface& interface() const { return *interface_.get(); }

  Controls& controls() { return controls_; }

  GLFWwindow* window() { return window_; }
  const GLFWwindow* window() const { return window_; }

  int width() const { return width_; }
  int height() const { return height_; }

  index_t fps() const { return fps_; }
  void set_fps( index_t fps ) { fps_ = fps; }

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

  std::shared_ptr<Interface> interface_;

  bool updated_;

  index_t fps_;
};

} // graphics

} // avro

#endif
