#ifndef avro_LIB_GRAPHICS_USER_INTERFACE_H_
#define avro_LIB_GRAPHICS_USER_INTERFACE_H_

#include "common/types.h"

#include "graphics/listener.h"

#include <imgui/imgui.h>

#include <memory>
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

protected:
  GLFW_Window& window_;
  const ImGuiIO& context_;
  Listener* listener_;
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

private:
  Visualizer& application_;

};

} // graphics

} // avro

#endif
