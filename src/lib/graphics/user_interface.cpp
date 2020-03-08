#include "common/error.h"

#include "graphics/user_interface.h"
#include "graphics/window.h"

#include <imgui/GL/imgui_impl_glfw.h>
#include <imgui/GL/imgui_impl_opengl3.h>

namespace avro
{

namespace graphics
{

Widget::Widget( const GLFW_Window& window ) :
  window_(window),
  context_(window.interface().context())
{}

Interface::Interface( GLFW_Window& window ) :
  window_(window)
{
  initialize();
}

void
Interface::initialize()
{
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  //context_ = ImGui::GetIO();

  // setup Dear ImGui style
  ImGui::StyleColorsDark();
  //ImGui::StyleColorsClassic();

  const char* glsl_version = "#version 410";

  // Setup Platform/Renderer bindings
  ImGui_ImplGlfw_InitForOpenGL(window_.window(), true);
  ImGui_ImplOpenGL3_Init(glsl_version);
}

void
Interface::begin_draw() const
{
  for (index_t k=0;k<widgets_.size();k++)
    widgets_[k]->begin_draw();
}

void
Interface::end_draw() const
{
  for (index_t k=0;k<widgets_.size();k++)
    widgets_[k]->end_draw();
}


Toolbar::Toolbar( const GLFW_Window& window ) :
  Widget(window)
{}

void
Toolbar::begin_draw() const
{
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();

  ImGuiWindowFlags window_flags = 0;
  window_flags |= ImGuiWindowFlags_NoCollapse;

  bool active = true;
  ImGui::Begin("Plot controls",&active,ImGuiWindowFlags_MenuBar);
  {
    if (ImGui::BeginMenuBar())
    {
      if (ImGui::BeginMenu("File"))
      {
          if (ImGui::MenuItem("Open..", "Ctrl+O")) { /* Do stuff */ }
          if (ImGui::MenuItem("Save", "Ctrl+S"))
          {
            printf("save file!\n");
          }
          if (ImGui::MenuItem("Close", "Ctrl+W"))  { active = false; }
          ImGui::EndMenu();
      }
      if (ImGui::BeginMenu("Edit"))
      {
          if (ImGui::MenuItem("Colors", NULL)) { /* Do stuff */ }
          ImGui::EndMenu();
      }
      if (ImGui::BeginMenu("View"))
      {
          if (ImGui::MenuItem("2D window", NULL)) { /* Do stuff */ }
          ImGui::EndMenu();
      }
      if (ImGui::BeginMenu("Help"))
      {
          if (ImGui::MenuItem("Usage", "Ctrl+H")) { /* Do stuff */ }
          ImGui::EndMenu();
      }
      ImGui::EndMenuBar();
    }

    if (ImGui::CollapsingHeader("Plots"))
    {
      if (ImGui::TreeNode("Plot 1"))
      {
        ImGui::Text("hello!");

        if (ImGui::TreeNode("Face 1"))
        {
          if (ImGui::TreeNode("Edge1"))
          {
            if (ImGui::TreeNode("Node 2"))
            {

              ImGui::TreePop();
            }
            ImGui::TreePop();
          }
          ImGui::TreePop();
        }
        ImGui::TreePop();
      }
      if (ImGui::TreeNode("Plot 2"))
      {
        ImGui::TreePop();
      }
    }

    if (ImGui::CollapsingHeader("Help"))
    {
      ImGui::Text("PROGRAMMER GUIDE:");
      ImGui::BulletText("Please see the ShowDemoWindow() code in imgui_demo.cpp. <- you are here!");
      ImGui::BulletText("Please see the comments in imgui.cpp.");
      ImGui::BulletText("Please see the examples/ application.");
      ImGui::BulletText("Enable 'io.ConfigFlags |= NavEnableKeyboard' for keyboard controls.");
      ImGui::BulletText("Enable 'io.ConfigFlags |= NavEnableGamepad' for gamepad controls.");
      ImGui::Separator();

      ImGui::Text("USER GUIDE:");
      //showUserGuide();
    }
  }

  ImGui::End();
}

void
Toolbar::end_draw() const
{
  ImGui::Render();
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

} // graphics

} // avro
