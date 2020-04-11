#include "../bin/programs.h"

#include "common/error.h"

#include "graphics/application.h"
#include "graphics/user_interface.h"
#include "graphics/window.h"

#include "library/factory.h"
#include "library/library.h"

#include "mesh/mesh.h"

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


Toolbar::Toolbar( GLFW_Window& window , Visualizer& app ) :
  Widget(window),
  application_(app)
{}

void
Toolbar::begin_draw()
{
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();

  // the global library that we will add to
  Library* lib = Library::get();

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

      ImGui::SameLine();
      if (ImGui::Button("Load"))
      {
        std::string ext = get_file_ext( items[file_current] );
        if (ext=="mesh")
        {
          lib->add_mesh( items[file_current] , listener_->pwd() );
        }
      }
    }

    if (ImGui::CollapsingHeader("Plots"))
    {
      std::vector<const char*> items;

      ImGui::SetNextItemWidth(100);
      static int mesh_current_plot = 0;
      const std::vector<std::string>& meshes = lib->meshes();
      items.resize( meshes.size() );
      for (index_t k=0;k<meshes.size();k++)
        items[k] = meshes[k].c_str();
      ImGui::Combo("Plot Mesh",&mesh_current_plot, items.data() , items.size() );

      ImGui::SameLine();
      if (ImGui::Button("Load Plot"))
      {
        // get the mesh
        std::string path = lib->meshname2file(meshes[mesh_current_plot]);
        if (path!="n/a")
        {
          // todo generalize this to other mesh types
          // i.e. determine what type of mesh this is first
          typedef Simplex type;

          // load the mesh
          std::shared_ptr<Topology<type>> ptopology = nullptr;
          std::shared_ptr<Mesh> pmesh = library::get_mesh<type>(path,ptopology);
          Topology<type>& topology = *ptopology.get();

          // load the topology and write to the visualizer
          // we first need someone to hold onto the mesh pointer (the library)
          lib->add_mesh_ptr(pmesh);
          application_.add_topology(topology);
          application_.restart() = true;
        }
      }

      index_t counter = 0;
      for (index_t k=0;k<window_.nb_scene();k++)
      {
        // skip this if we need to restart
        if (application_.restart()) break;

        const SceneGraph& scene = window_.scene(k);
        const json& menu = scene.menu();

        const std::vector<std::string>& primitives = menu["primitives"];
        for (index_t j=0;j<primitives.size();j++)
        {
          json plot = json::parse(primitives[j]);

          std::string label = "Plot " + std::to_string(counter++);
          if (ImGui::TreeNode(label.c_str()))
          {
            if (!window_.scene(k).primitive(j).hidden())
            {
              std::string hide_label = "Hide" + std::to_string(counter-1);
              if (ImGui::Button(hide_label.c_str()))
              {
                for (index_t i=0;i<window_.scene(j).nb_primitives();i++)
                {
                  window_.scene(k).primitive(j).hide();
                }
              }
            }
            else
            {
              std::string show_label = "Show" + std::to_string(counter-1);
              if (ImGui::Button(show_label.c_str()))
              {
                for (index_t i=0;i<window_.scene(j).nb_primitives();i++)
                {
                  window_.scene(k).primitive(j).show();
                }
              }
            }

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

/*                  if (hide)
                  {
                    primitive->triangles_on() = false;
                    primitive->edges_on() = false;
                  }*/

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

    if (ImGui::CollapsingHeader("Adapt"))
    {
      std::vector<const char*> items;

      ImGui::SetNextItemWidth(100);
      static int mesh_current = 0;
      const std::vector<std::string>& meshes = lib->meshes();
      items.resize( meshes.size() );
      for (index_t k=0;k<meshes.size();k++)
        items[k] = meshes[k].c_str();
      ImGui::Combo("Mesh",&mesh_current, items.data() , items.size() );

      ImGui::SameLine();
      ImGui::SetNextItemWidth(5*width);
      static char iter_str[5] = "1";
      ImGui::InputText("Iter.", iter_str , IM_ARRAYSIZE(iter_str));
      int nb_adapt = atoi(iter_str);

      ImGui::SameLine();
      static char output_name[128] = "[output].mesh";
      ImGui::SetNextItemWidth(100);
      ImGui::InputText("Output:", output_name, IM_ARRAYSIZE(output_name));

      ImGui::SetNextItemWidth(100);
      static int geometry_current = 0;
      const std::vector<std::string>& geometries = lib->geometries();
      items.resize( geometries.size() );
      for (index_t k=0;k<geometries.size();k++)
        items[k] = geometries[k].c_str();
      ImGui::Combo("Geometry",&geometry_current, items.data() , items.size() );

      coord_t dim = 3;

      ImGui::SetNextItemWidth(100);
      static int metric_current = 0;
      const std::vector<std::string>& metrics = lib->metrics();
      items.resize( metrics.size() );
      for (index_t k=0;k<metrics.size();k++)
        items[k] = metrics[k].c_str();
      ImGui::Combo("Metric",&metric_current, items.data() , items.size() );

      std::string mesh_name = meshes[mesh_current];
      if (ImGui::Button("Compute"))
      {
        ImGui::OpenPopup("adapt-mesh");

      }
      if (ImGui::BeginPopupModal("adapt-mesh", NULL, ImGuiWindowFlags_MenuBar))
      {
        if (mesh_name=="CKF")
        {
          for (index_t d=0;d<dim;d++)
            mesh_name += "-3";
        }
        ImGui::Text("Adapt mesh %s, metric %s, geometry %s?\n",mesh_name.c_str(),metrics[metric_current].c_str(),geometries[geometry_current].c_str());
        if (ImGui::Button("Adapt"))
        {

          ImGui::CloseCurrentPopup();
          std::string iter_cmd = "nb_iter="+std::to_string(nb_adapt);
          const char* command[] = {mesh_name.c_str(),geometries[geometry_current].c_str(),metrics[metric_current].c_str(),output_name,iter_cmd.c_str()};
          int result = programs::adapt(5,command);
        }
        ImGui::EndPopup();
      }
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

    if (ImGui::CollapsingHeader("Tools"))
    {
      ImGui::SetNextItemWidth(100);

      static int fps = 60;
      ImGui::SliderInt("fps",&fps,5,120,"%3d");
      window_.set_fps( fps );
    }

    if (ImGui::CollapsingHeader("Help"))
    {
      ImGui::Text("avro (c) Philip Caplan 2019-2020\ncaplan.philip@gmail.com");
      ImGui::BulletText("rotatee: hold ctrl (cmd) and click/move mouse");
      ImGui::BulletText("zoom: hold shift key and click/move mouse");
      ImGui::BulletText("pan: click/move mouse");
      ImGui::Separator();
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
