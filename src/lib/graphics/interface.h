#ifndef URSA_LIB_GRAPHICS_INTERFACE_H_
#define URSA_LIB_GRAPHICS_INTERFACE_H_

#include <imgui/imgui.h>

namespace ursa
{

namespace graphics
{

class Window;

class InterfaceManager
{
public:
  InterfaceManager();

  void initialize();

  const ImGuiIO& context() const { return context_; }
  ImGuiIO& context() { return context_; }

private:
  ImGuiIO context_;
};

class Interface
{
public:
  Interface( Window& window );
  virtual ~Interface() {}

  virtual void show() = 0;
  void render();

private:
  InterfaceManager& manager_;
  Window& window_;
};

class BasicInterface : public Interface
{
public:
  BasicInterface( Window& window );

  void show(); // file, edit, view, etc...
};

class PlotTree : public Interface
{
public:
  PlotTree( Window& window );
  void show();
};

};

};

#endif
