#ifndef URSA_LIB_GRAPHICS_WINDOW_H_
#define URSA_LIB_GRAPHICS_WINDOW_H_

#include "graphics/gl.h"
#include "graphics/math.h"
#include "graphics/scene.h"

#include <memory>
#include <string>

namespace ursa
{

namespace graphics
{

class Plot;
class Plotter;
class Scene;

class Window
{
public:
  typedef std::shared_ptr<Plot> Plot_ptr;

  Window( const std::string& title , Plotter* plotter );
  virtual ~Window();

  virtual void draw();
  virtual void write();

  void run();

  std::string title() const { return title_; }

  Scene& scene() { return scene_; }
  const Scene& scene() const { return scene_; }

  void attach( Plot_ptr plot );

  GLFWwindow* window() { return window_; }

  Plotter* plotter() { return plotter_; }
  const Plotter* plotter() const { return plotter_; }

private:
  std::string title_;
  Scene scene_;
  Plotter* plotter_;
  GLFWwindow* window_;

  int width_ = 650;
  int height_ = 400;

  std::vector<Plot_ptr> plot_;

  mat3 mvp_;
};

} // graphics

} // ursa

#endif
