#ifndef avro_LIB_GRAPHICS_USER_INTERFACE_H_
#define avro_LIB_GRAPHICS_USER_INTERFACE_H_

#include "common/types.h"

#include <imgui/imgui.h>

#include <memory>
#include <vector>

namespace avro
{

namespace graphics
{

class GLFW_Window;

class Widget
{
protected:
  Widget( const GLFW_Window& window );
  virtual ~Widget() {}

public:
  virtual void begin_draw() const = 0;
  virtual void end_draw() const = 0;

protected:
  const GLFW_Window& window_;
  const ImGuiIO& context_;
};

class Interface
{
public:
  Interface( GLFW_Window& window );

  void initialize();

  void begin_draw() const;
  void end_draw() const;

  void add_widget( std::shared_ptr<Widget> widget )
  {
    widgets_.push_back(widget);
  }

  const ImGuiIO& context() const { return context_; }
  ImGuiIO& context() { return context_; }

private:
  GLFW_Window& window_;
  ImGuiIO context_;

  std::vector<std::shared_ptr<Widget>> widgets_;
};

class Toolbar : public Widget
{
public:
  Toolbar( const GLFW_Window& window );

  void begin_draw() const;
  void end_draw() const;

private:
  // nothing yet
};

} // graphics

} // avro

#endif