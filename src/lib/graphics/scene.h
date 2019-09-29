#ifndef URSA_LIB_GRAPHICS_SCENE_H_
#define URSA_LIB_GRAPHICS_SCENE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace ursa
{

namespace graphics
{

class Plot;
class Shader;

class Scene
{
  typedef std::shared_ptr<Plot> Plot_ptr;
  typedef std::shared_ptr<Shader> Shader_ptr;

public:
  Scene();

  void render();

private:
  std::map< std::string , Plot_ptr > plot_;

  std::vector< Shader_ptr > shaders_;
};

} // graphics

} // ursa

#endif
