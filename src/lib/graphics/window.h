#ifndef URSA_LIB_GRAPHICS_WINDOW_H_
#define URSA_LIB_GRAPHICS_WINDOW_H_

#include "graphics/scene.h"

#include <string>

namespace ursa
{

namespace graphics
{

class Scene;

class Window
{
public:
  Window( const std::string& title );

  void initialize();

  std::string title() const { return title_; }

  Scene& scene() { return scene_; }
  const Scene& scene() const { return scene_; }

private:
  std::string title_;
  Scene scene_;
};

} // graphics

} // ursa

#endif
