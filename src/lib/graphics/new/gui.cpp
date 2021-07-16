#include "graphics/new/gui.h"
#include "graphics/new/window.h"

#include "graphics/gl.h"

#include <imgui/GL/imgui_impl_glfw.h>
#include <imgui/GL/imgui_impl_opengl3.h>

namespace avro
{

namespace graphics
{

GUI::GUI( Window& window ) :
 window_(window),
 count_(5)
{
  window_.set_gui(this);

  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  //context_ = ImGui::GetIO();

  // setup Dear ImGui style
  ImGui::StyleColorsDark();
  //ImGui::StyleColorsClassic();

#if AVRO_HEADLESS_GRAPHICS
  const char* glsl_version = "#version 330";
#else
  const char* glsl_version = "#version 410";
  //const char* glsl_version = "#version 330"; // again, I would like this to be 4.1
#endif

  // setup Platform/Renderer bindings
  ImGui_ImplGlfw_InitForOpenGL(window_.window(), true);
  ImGui_ImplOpenGL3_Init(glsl_version);
}

void
GUI::begin_draw() {
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();
}

void
GUI::draw() {

  begin_draw();

  //bool capture_mouse = ImGui::GetIO().WantCaptureMouse;
  ImGuiWindowFlags window_flags = 0;
  //window_flags |= ImGuiWindowFlags_NoCollapse;
  window_flags |= ImGuiWindowFlags_NoResize;
  //window_flags |= ImGuiWindowFlags_AlwaysHorizontalScrollbar;

  bool capture_mouse = ImGui::GetIO().WantCaptureMouse;
  if (capture_mouse) window_.needs_drawing(true);

  // set the controls at the top-left of the window
  ImGui::SetNextWindowPos( ImVec2( 0 , 0 ) );
  ImGui::SetNextWindowSize( ImVec2( 200 , window_.height() ) );

  std::vector<std::string> entities = {"Faces","Edges","Nodes"};

  bool active = true;
  ImGui::SetNextItemWidth(200);

  if (ImGui::Begin("Controls",&active,window_flags))
  {
    ImGui::Text("hello");

    index_t nb_plots = 3;
    for (index_t i = 0; i < nb_plots; i++) {

      std::string label = "Plot" + std::to_string(i);
      if (ImGui::CollapsingHeader(label.c_str())) {


      }
    }

    ImGui::Text("FPS %.1f", ImGui::GetIO().Framerate);

    ImGui::End();
  }

  ImGui::Render();
  window_.draw();
  #if AVRO_HEADLESS_GRAPHICS == 0
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
  #endif

  glfwSwapBuffers(window_.window());
}


} // graphics

} // avro
