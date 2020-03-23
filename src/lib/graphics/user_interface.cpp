#include "common/error.h"

#include "graphics/user_interface.h"
#include "graphics/window.h"

#include <imgui/GL/imgui_impl_glfw.h>
#include <imgui/GL/imgui_impl_opengl3.h>

namespace avro
{

namespace graphics
{

Widget::Widget( GLFW_Window& window ) :
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
Interface::begin_draw()
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


Toolbar::Toolbar( GLFW_Window& window ) :
  Widget(window)
{}

void
Toolbar::begin_draw()
{
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();

  // determine if imgui wants the mouse or if we should send it to the trackball
  bool capture_mouse = ImGui::GetIO().WantCaptureMouse;
  if (capture_mouse)
    window_.controls().disable();
  else
    window_.controls().enable();

  ImGuiWindowFlags window_flags = 0;
  window_flags |= ImGuiWindowFlags_NoCollapse;

  std::vector<std::string> entities = {"Volumes","Faces","Edges","Nodes"};

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
      for (index_t k=0;k<window_.nb_scene();k++)
      {
        const SceneGraph& scene = window_.scene(k);
        const json& menu = scene.menu();

        const std::vector<std::string>& primitives = menu["primitives"];
        for (index_t j=0;j<primitives.size();j++)
        {
          const json& plot = json::parse(primitives[k]);

          std::string label = "Plot " + std::to_string(j);
          if (ImGui::TreeNode(label.c_str()))
          {
            for (index_t i=0;i<entities.size();i++)
            {
              if (ImGui::TreeNode(entities[i].c_str()))
              {
                std::vector<std::string> entity = plot[entities[i]];
                for (index_t m=0;m<entity.size();m++)
                {
                  static bool viz = true;
                  static float alpha = 1.0f;
                  static bool edg = true;
                  ImGui::Text(entity[m].c_str());
                  ImGui::SameLine();
                  ImGui::Checkbox("viz",&viz);
                  ImGui::SameLine();
                  ImGui::Checkbox("edg",&edg);
                  ImGui::SameLine();
                  ImGui::PushItemWidth(50);
                  ImGui::SliderFloat("alpha",&alpha,0.1f,1.0f,"%.2f");

                  unsigned long address = std::stoul(entity[m],0,16);
                  Primitive* primitive = (Primitive*) address;

                  primitive->triangles_on() = viz;
                  primitive->edges_on() = edg;
                  primitive->set_transparency(alpha);
                }
                ImGui::TreePop();
              }
            }
            ImGui::TreePop();
          }
        }

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
