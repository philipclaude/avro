#include "common/tools.h"

#include "graphics/interface.h"
#include "graphics/plotter.h"
#include "graphics/window.h"

#include <imgui/GL/imgui_impl_glfw.h>
#include <imgui/GL/imgui_impl_opengl3.h>

namespace avro
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

// Helper to display basic user controls.
void showUserGuide()
{
    ImGuiIO& io = ImGui::GetIO();
    ImGui::BulletText("Double-click on title bar to collapse window.");
    ImGui::BulletText("Click and drag on lower corner to resize window\n(double-click to auto fit window to its contents).");
    if (io.ConfigWindowsMoveFromTitleBarOnly)
        ImGui::BulletText("Click and drag on title bar to move window.");
    else
        ImGui::BulletText("Click and drag on any empty space to move window.");
    ImGui::BulletText("TAB/SHIFT+TAB to cycle through keyboard editable fields.");
    ImGui::BulletText("CTRL+Click on a slider or drag box to input value as text.");
    if (io.FontAllowUserScaling)
        ImGui::BulletText("CTRL+Mouse Wheel to zoom window contents.");
    ImGui::BulletText("Mouse Wheel to scroll.");
    ImGui::BulletText("While editing text:\n");
    ImGui::Indent();
    ImGui::BulletText("Hold SHIFT or use mouse to select text.");
    ImGui::BulletText("CTRL+Left/Right to word jump.");
    ImGui::BulletText("CTRL+A or double-click to select all.");
    ImGui::BulletText("CTRL+X,CTRL+C,CTRL+V to use clipboard.");
    ImGui::BulletText("CTRL+Z,CTRL+Y to undo/redo.");
    ImGui::BulletText("ESCAPE to revert.");
    ImGui::BulletText("You can apply arithmetic operators +,*,/ on numerical values.\nUse +- to subtract.");
    ImGui::Unindent();
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
            window_.save("test.eps");
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
      showUserGuide();
    }
  }

  ImGui::End();

  ImGui::Render();
}

} // graphics

} // avro
