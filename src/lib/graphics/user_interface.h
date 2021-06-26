//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_GRAPHICS_USER_INTERFACE_H_
#define avro_LIB_GRAPHICS_USER_INTERFACE_H_

#include "types.h"

#include "graphics/listener.h"

#include <imgui/imgui.h>

#include <memory>
#include <set>
#include <vector>

namespace avro
{

namespace graphics
{

class GLFW_Window;
class Visualizer;

class Widget
{
protected:
  Widget( GLFW_Window& window );
  virtual ~Widget() {}

public:
  virtual void begin_draw() = 0;
  virtual void end_draw() const = 0;

  void set_listener( Listener* listener )
  { listener_ = listener; }

  bool& active() { return active_; }

protected:
  GLFW_Window& window_;
  const ImGuiIO& context_;
  Listener* listener_;
  bool active_;
};

class Interface
{
public:
  Interface( GLFW_Window& window , Listener& listener );

  void initialize();

  void begin_draw();
  void end_draw() const;

  void add_widget( std::shared_ptr<Widget> widget )
  {
    widgets_.push_back(widget);
    widget->set_listener(&listener_);
  }

  const ImGuiIO& context() const { return context_; }
  ImGuiIO& context() { return context_; }

  bool& active() { return widgets_[0]->active(); }

private:
  GLFW_Window& window_;
  ImGuiIO context_;

  std::vector<std::shared_ptr<Widget>> widgets_;

  Listener& listener_;
};

class Toolbar : public Widget
{
public:
  Toolbar( GLFW_Window& window , Visualizer& app );

  void begin_draw();
  void end_draw() const;

  void save_eps(bool& open) const;

private:
  Visualizer& application_;

  std::map<std::string,std::pair<index_t,index_t>> primitives_;

};

} // graphics

} // avro

#endif
