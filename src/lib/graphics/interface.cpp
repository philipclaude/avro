#include "common/tools.h"

#include "graphics/interface.h"
#include "graphics/plotter.h"
#include "graphics/window.h"

#include <imgui/GL/imgui_impl_glfw.h>
#include <imgui/GL/imgui_impl_opengl3.h>

namespace ursa
{

namespace graphics
{

InterfaceManager::InterfaceManager()
{
  initialize();
}

void
InterfaceManager::initialize()
{
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  context_ = ImGui::GetIO();

  // Setup Dear ImGui style
  ImGui::StyleColorsDark();
  //ImGui::StyleColorsClassic();
}

Interface::Interface( Window& window ) :
  manager_(window.plotter()->manager()),
  window_(window)
{
  UNUSED(manager_);
}

void
Interface::render()
{
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

BasicInterface::BasicInterface( Window& window ) :
  Interface(window)
{}

void
BasicInterface::show()
{
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();

/*
  bool show_demo_window = true;
  ImGui::ShowDemoWindow(&show_demo_window);
*/

  // initialize the ImGui render before we actually render with the gl
  ImGui::Render();
}

PlotTree::PlotTree( Window& window ) :
  Interface(window)
{}

void
PlotTree::show()
{
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();

  ImGuiWindowFlags window_flags = 0;
  window_flags |= ImGuiWindowFlags_NoCollapse;

  bool p_open = false;
  if (!ImGui::Begin("primitive tree",&p_open,window_flags))
  {
      // Early out if the window is collapsed, as an optimization.
      ImGui::End();
      return;
  }

  {
      if (ImGui::BeginMenu("File"))
      {
        bool show_file = true;
        bool show_edit = true;
        bool show_view = true;
        bool show_help = true;
        ImGui::MenuItem("File",NULL,&show_file);
        ImGui::MenuItem("Edit",NULL,&show_edit);
        ImGui::MenuItem("View",NULL,&show_view);
        ImGui::MenuItem("Help",NULL,&show_help);

        ImGui::EndMenu();
      }
      ImGui::SameLine();
      if (ImGui::BeginMenu("Edit"))
      {
        bool show_file = true;
        bool show_edit = true;
        bool show_view = true;
        bool show_help = true;
        ImGui::MenuItem("File",NULL,&show_file);
        ImGui::MenuItem("Edit",NULL,&show_edit);
        ImGui::MenuItem("View",NULL,&show_view);
        ImGui::MenuItem("Help",NULL,&show_help);

        ImGui::EndMenu();
      }

      if (ImGui::TreeNode("Style"))
      {
          ImGui::ShowStyleEditor();
          ImGui::TreePop();
          ImGui::Separator();
      }
  }

  ImGui::End();

  ImGui::Render();
}

} // graphics

} // ursa
