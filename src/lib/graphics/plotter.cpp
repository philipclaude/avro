#include "common/error.h"

#include <glad/glad.h>

#include "graphics/gl.h"
#include "graphics/plotter.h"
#include "graphics/shader.h"
#include "graphics/window.h"

#include <imgui/GL/imgui_impl_glfw.h>
#include <imgui/GL/imgui_impl_opengl3.h>

namespace luna
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
  luna_assert_msg( glfwInit() , "problem initializing OpenGL!" );

  // set the version
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  // create and make the window current
  createWindow( "luna plot" );
  GLFWwindow* window = window_["luna plot"]->window();
  glfwMakeContextCurrent(window);

  // load gl stuff and print info
  gladLoadGL();
  dumpGLInfo();

  // create all the shaders
  shader_.insert( {"basic" , std::make_shared<ShaderProgram>( "basic" ) } );
  shader_.insert( {"edge" , std::make_shared<ShaderProgram>( "edge" ) } );

  shader_.insert( {"wv" , std::make_shared<ShaderProgram>( "wv" ) } );
  // TODO the rest of the shaders...

  const char* glsl_version = "#version 410";
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO& io = ImGui::GetIO(); (void)io;
  //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
  //io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls

  // Setup Dear ImGui style
  ImGui::StyleColorsDark();
  //ImGui::StyleColorsClassic();

  // Setup Platform/Renderer bindings
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL3_Init(glsl_version);

  bool show_demo_window = true;
  bool show_another_window = false;

/*
  ImGui::Begin("Another Window", &show_another_window);   // Pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
  ImGui::Text("Hello from another window!");
  if (ImGui::Button("Close Me"))
      show_another_window = false;
  ImGui::End();
*/

}

void
Plotter::run()
{
  for (std::map<std::string,Window_ptr>::iterator it=window_.begin();it!=window_.end();it++)
    it->second->run();
}

Plotter::~Plotter()
{
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();

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
  if (name=="main") return window("luna plot");

  luna_assert( window_.find(name)!=window_.end() );
  return *window_[name].get();
}

ShaderProgram&
Plotter::shader( const std::string& name )
{
  luna_assert( shader_.find(name)!=shader_.end() );
  return *shader_.at(name).get();
}

} // graphics

} // luna
