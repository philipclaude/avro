#ifndef URSA_LIB_GRAPHICS_PLOTTER_H_
#define URSA_LIB_GRAPHICS_PLOTTER_H_

#include "graphics/client.h"
#include "graphics/server.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace ursa
{

namespace graphics
{

class Window;
class ShaderProgram;

class Plotter
{
  typedef std::shared_ptr<Window> Window_ptr;
  typedef std::shared_ptr<ShaderProgram> Shader_ptr;

public:

  Plotter();
  ~Plotter();

  void initialize();

  void run();

  void createWindow( const std::string& title );

  Window& window( const std::string& name );

  ShaderProgram& shader( const std::string& name );
  const ShaderProgram& shader( const std::string& name ) const;


private:
  std::map<std::string,Window_ptr> window_;
  std::map<std::string,Shader_ptr> shader_;
};

} // graphics

} // ursa

#endif
