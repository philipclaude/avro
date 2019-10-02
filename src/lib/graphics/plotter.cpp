#include "common/error.h"

#include <glad/glad.h>

#include "graphics/gl.h"
#include "graphics/plotter.h"
#include "graphics/shader.h"
#include "graphics/window.h"

namespace ursa
{

namespace graphics
{

Plotter::Plotter()
{
  initialize();
}

void
Plotter::initialize()
{
  // initialize OpenGL
  ursa_assert_msg( glfwInit() , "problem initializing OpenGL!" );

  // set the version
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  // create and make the window current
  createWindow( "ursa plotter" );
  GLFWwindow* window = window_["ursa plotter"]->window();
  glfwMakeContextCurrent(window);

  // load gl stuff and print info
  gladLoadGL();
  dumpGLInfo();

  // create all the shaders
  shader_.insert( {"basic" , std::make_shared<ShaderProgram>( "basic" ) } );
  shader_.insert( {"edge" , std::make_shared<ShaderProgram>( "edge" ) } );
  // TODO the rest of the shaders...
}

void
Plotter::run()
{
  for (std::map<std::string,Window_ptr>::iterator it=window_.begin();it!=window_.end();it++)
    it->second->run();
}

Plotter::~Plotter()
{
  glfwTerminate();
}

void
Plotter::createWindow( const std::string& title )
{
  std::shared_ptr<Window> w = std::make_shared<Window>(title,this);
  window_.insert( {title  , w } );
}

Window&
Plotter::window( const std::string& name )
{
  if (name=="main") return window("ursa plotter");

  ursa_assert( window_.find(name)!=window_.end() );
  return *window_[name].get();
}

ShaderProgram&
Plotter::shader( const std::string& name )
{
  ursa_assert( shader_.find(name)!=shader_.end() );
  return *shader_.at(name).get();
}

} // graphics

} // ursa
