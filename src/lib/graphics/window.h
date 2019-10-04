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

  void setMatrices();

  virtual void draw();
  virtual void write();

  void run();

  std::string title() const { return title_; }

  Scene& scene() { return scene_; }
  const Scene& scene() const { return scene_; }

  void attach( Plot_ptr plot );

  GLFWwindow* window() { return window_; }

  const mat4& mvp() const { return mvp_; }
  mat4& mvp() { return mvp_; }

  const mat4& viewMatrix() const { return viewMatrix_; }
  mat4& viewMatrix() { return viewMatrix_; }

  const mat4& projMatrix() const { return projMatrix_; }
  mat4& projMatrix() { return projMatrix_; }

  const mat4& modelMatrix() const { return modelMatrix_; }
  mat4& modelMatrix() { return modelMatrix_; }

  Plotter* plotter() { return plotter_; }
  const Plotter* plotter() const { return plotter_; }

private:
  std::string title_;
  Scene scene_;
  Plotter* plotter_;
  GLFWwindow* window_;

  int width_ = 1024;//650;
  int height_ = 640;//400;

  std::vector<Plot_ptr> plot_;

  mat4 mvp_;
  mat4 viewMatrix_;
  mat4 projMatrix_;
  mat4 modelMatrix_;

  float fov_ = 45.0f;
  float speed_ = 3.0f;
  float mouseSpeed_ = 0.005f;

  vec3 position_;

  float angles_[2];
};

} // graphics

} // ursa

#endif
