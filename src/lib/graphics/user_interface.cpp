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
  context_(window.interface().context()),
  listener_(nullptr)
{}

Interface::Interface( GLFW_Window& window , Listener& listener ) :
  window_(window),
  listener_(listener)
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
  ImGui::SetNextItemWidth(200);
  ImGui::Begin("Controls",&active,ImGuiWindowFlags_MenuBar);
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

    ImGuiIO& io = ImGui::GetIO();
    ImFont* pfont = io.Fonts->Fonts[0];
    int width = pfont->FontSize;

    ImGui::SetNextItemWidth(100);
    if (ImGui::CollapsingHeader("I/O"))
    {
      json request,response;
      request["command"] = "ls";
      listener_->send(request,response);
      const std::vector<json>& directory = response["ls-response"];

      std::vector<std::string> dirs;
      std::vector<std::string> files;
      index_t lmax = 10;
      for (index_t k=0;k<directory.size();k++)
      {
        const std::string& s = directory[k]["entry"];
        if (directory[k]["type"] == "file")
          files.push_back( s );
        else if (directory[k]["type"] == "dir")
          dirs.push_back( s );
        lmax = std::max ( lmax , s.size() );
      }

      std::vector<const char*> items;
      items.push_back("[pwd]");
      for (index_t k=0;k<dirs.size();k++)
        items.push_back( dirs[k].c_str() );



      static int dir_current = 0;
      ImGui::SetNextItemWidth(lmax*width);
      ImGui::Combo("Directory", &dir_current, items.data() , items.size() );

      if (dir_current!=0)
      {
        request["command"] = "cd";
        request["data"]    = dirs[dir_current-1];
        listener_->send(request,response);
        dir_current = 0;
      }

      items.clear();
      for (index_t k=0;k<files.size();k++)
        items.push_back( files[k].c_str() );

      static int file_current = 0;
      ImGui::SetNextItemWidth(lmax*width);
      ImGui::Combo("File",&file_current, items.data() , items.size() );

      static char str0[128] = "[id name]";
      ImGui::SetNextItemWidth(lmax*width);
      ImGui::InputText("", str0, IM_ARRAYSIZE(str0));
      ImGui::SameLine();
      ImGui::Button("Load");
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
                  static float alpha = 1.0f;

                  unsigned long address = std::stoul(entity[m],0,16);
                  Primitive* primitive = (Primitive*) address;

                  ImGui::Text("%s",entity[m].c_str());
                  ImGui::SameLine();
                  std::string viz_label = "viz-" + std::to_string(m);
                  ImGui::Checkbox(viz_label.c_str(), &primitive->triangles_on() );
                  ImGui::SameLine();
                  std::string edg_label = "edg-" + std::to_string(m);
                  ImGui::Checkbox(edg_label.c_str(), &primitive->edges_on() );
                  ImGui::SameLine();
                  ImGui::PushItemWidth(50);
                  std::string alpha_label = "alpha-" + std::to_string(m);
                  ImGui::SliderFloat(alpha_label.c_str(),&primitive->transparency(),0.1f,1.0f,"%.2f");

                }
                ImGui::TreePop();
              }
            }
            ImGui::TreePop();
          }
        }

      }
    }

    if (ImGui::CollapsingHeader("Tools"))
    {
      ImGui::SetNextItemWidth(100);

      static int fps = 60;
      ImGui::SliderInt("fps",&fps,5,120,"%3d");
      window_.set_fps( fps );
    }

    if (ImGui::CollapsingHeader("Adapt"))
    {
      ImGui::SetNextItemWidth(100);
      static int mesh_current = 0;
      const char* meshes[] = {"mesh0","mesh1","mesh2"};
      ImGui::Combo("Mesh",&mesh_current, meshes , 3 );

      ImGui::SameLine();
      ImGui::SetNextItemWidth(5*width);
      static char iter_str[5] = "1";
      ImGui::InputText("Iter.", iter_str , IM_ARRAYSIZE(iter_str));
      int nb_adapt = atoi(iter_str);

      ImGui::SetNextItemWidth(100);
      static int geometry_current = 0;
      const char* geometries[] = {"geometry0","geometry1","geometry2"};
      ImGui::Combo("Geometry",&geometry_current, geometries , 3 );

      ImGui::SetNextItemWidth(100);
      static int metric_current = 0;
      const char* metrics[] = {"metric0","metric1","metric2"};
      ImGui::Combo("Metric",&metric_current, metrics , 3 );

      ImGui::Button("Load");
      ImGui::SameLine();
      ImGui::Button("Compute");
    }

    if (ImGui::CollapsingHeader("Voronoi"))
    {
      ImGui::SetNextItemWidth(100);
      static int mesh_current = 0;
      const char* meshes[] = {"mesh0","mesh1","mesh2"};
      ImGui::Combo("Mesh",&mesh_current, meshes , 3 );

      ImGui::SetNextItemWidth(100);
      static int sites_current = 0;
      const char* sites[] = {"vertices","centroids","random"};
      ImGui::Combo("Sites",&sites_current, sites , 3 );

      ImGui::Button("Load");
      ImGui::SameLine();
      ImGui::Button("Compute");
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
