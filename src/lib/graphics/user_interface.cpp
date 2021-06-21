//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "../bin/programs.h"

#include "common/error.h"
#include "common/tools.h"

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

#ifdef AVRO_WITH_GL

Widget::Widget( GLFW_Window& window ) :
  window_(window),
  context_(window.interface().context()),
  listener_(nullptr),
  active_(false)
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

#ifdef AVRO_HEADLESS_GRAPHICS
  const char* glsl_version = "#version 330";
#else
  //const char* glsl_version = "#version 410";
  const char* glsl_version = "#version 330"; // again, I would like this to be 4.1
#endif

  // setup Platform/Renderer bindings
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
parse_parameters( const std::string& builtin , const std::string& params )
{
  Library* lib = Library::get();

  std::cout << builtin << std::endl;
  std::vector<std::string> s = split(builtin,":");
  if (s[0]=="Mesh")
  {
    avro_assert( s.size() == 2 );
    s[1].erase(remove_if(s[1].begin(), s[1].end(), isspace), s[1].end());
    std::cout << s[1] << std::endl;
    std::shared_ptr<TopologyBase> ptopology;
    std::shared_ptr<Mesh> pmesh = library::get_mesh(s[1]+"-"+params,ptopology);
    char mesh_label[128];
    sprintf(mesh_label,"%p",(void*)pmesh.get());
    lib->add_mesh(mesh_label);
    lib->add_mesh_ptr(pmesh);
  }
  if (s[1]=="Metric")
  {
    avro_implement;
  }
}

void
Toolbar::save_eps(bool& open) const
{
  ImGui::SetNextWindowSize(ImVec2(250, 100));
  ImGui::Begin("Save file", &open);

  static char output_name[128] = ".eps";
  ImGui::SetNextItemWidth(150);
  ImGui::InputText("Filename", output_name, IM_ARRAYSIZE(output_name));

  if (ImGui::Button("Save"))
  {

    std::string filename(output_name);
    application_.main_window().save_eps( filename );

    open = false;
  }
  ImGui::SameLine();
  if (ImGui::Button("Cancel"))
  {
    open = false;
  }

  ImGui::End();

}

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
  {
    window_.controls().disable();
    window_.clip_controls().disable();
    active_ = true;
  }
  else
  {
    window_.controls().enable();
    window_.clip_controls().enable();
    active_ = false;
  }

  ImGuiWindowFlags window_flags = 0;
  window_flags |= ImGuiWindowFlags_NoCollapse;

  std::vector<std::string> entities = {"Volumes","Faces","Edges","Nodes"};

  static bool show_save_eps = false;

  if (show_save_eps) save_eps(show_save_eps);

  bool active = true;
  ImGui::SetNextItemWidth(200);
  ImGui::Begin("Controls",&active,ImGuiWindowFlags_MenuBar);
  {
    if (ImGui::BeginMenuBar())
    {
      if (ImGui::BeginMenu("File"))
      {
          if (ImGui::MenuItem("Open..", "Ctrl+O")) { /* Do stuff */ }
          if (ImGui::BeginMenu("Save"))
          {
            if (ImGui::MenuItem("EPS",NULL,&show_save_eps))
            {
              ImGui::Text("enter file name");
              ImGui::SameLine();
              ImGui::Button("save file!");
              ImGui::OpenPopup("save-window");
            }

            ImGui::MenuItem("PNG");
            ImGui::EndMenu();
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
      ImGui::EndMenuBar();
    }

    ImGuiIO& io = ImGui::GetIO();
    ImFont* pfont = io.Fonts->Fonts[0];
    int width = pfont->FontSize;

    ImGui::SetNextItemWidth(100);
    if (ImGui::CollapsingHeader("Library"))
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
      if (ImGui::Button("Add"))
      {
        std::string ext = get_file_ext( items[file_current] );
        if (ext=="mesh")
        {
          lib->add_mesh( items[file_current] , listener_->pwd() );
        }
        if (ext=="egads")
        {
          lib->add_geometry( items[file_current] , listener_->pwd() );
        }
      }

      // custom parameter loader
      ImGui::SetNextItemWidth(100);
      items.clear();
      items.push_back("Mesh: CKF-simplex");
      items.push_back("Mesh: CKF-cube");
      items.push_back("Metric: Linear");
      static int builtin_current = 0;
      ImGui::Combo("Built-in",&builtin_current,items.data(),items.size());

      ImGui::SameLine();
      static char param_name[128] = "params";
      ImGui::SetNextItemWidth(100);
      ImGui::InputText("Parameters", param_name, IM_ARRAYSIZE(param_name));
      ImGui::SameLine();
      if (ImGui::Button("Load"))
      {
        std::string builtin(items[builtin_current]);
        std::string params(param_name);
        parse_parameters( builtin , params );
        printf("parse parameters and load mesh/metric/geometry!\n");
      }

    }

    if (ImGui::CollapsingHeader("Plots"))
    {
      primitives_.clear();
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
        std::string mesh_name;
        if (path!="n/a")
          mesh_name = path;
        else
          mesh_name = meshes[mesh_current_plot];

        // load the mesh
        std::shared_ptr<Mesh> pmesh = nullptr;
        std::shared_ptr<TopologyBase> ptopology = nullptr;
        pmesh = library::get_mesh(mesh_name,ptopology);
        TopologyBase& topology = *ptopology.get();
        application_.add_topology(topology);

        // load the topology and write to the visualizer
        // we first need someone to hold onto the mesh pointer (the library)
        lib->add_mesh_ptr(pmesh);
        application_.restart() = true;
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

            // fields dropdown
            ImGui::SetNextItemWidth(100);
            static int current_field = 0;
            std::vector<std::string> field_names = plot["fields"];
            std::vector<std::string> field_ids = plot["fields_id"];
            avro_assert( field_names.size() == field_ids.size() );
            items.resize( field_names.size() + 1 );
            items[0] = "none";
            for (index_t k=0;k<field_names.size();k++)
              items[k+1] = field_names[k].c_str();

            real_t ulim[2] = {0,1};
            ImGui::Combo("Fields",&current_field,items.data(),items.size());
            ImGui::SameLine();
            if (ImGui::Button("Load"))
            {
              std::string label(items[current_field]);

              if (current_field==0)
              {
                // reset to no colors
                window_.scene(k).primitive(j).set_active( "0x0" , 0 );
              }
              else
              {
                std::vector<std::string> s = split(field_ids[current_field-1],"-");
                avro_assert( s.size()==2 );
                index_t rank = atoi(s[1].c_str());

                const Fields& fields = window_.scene(k).primitive(j).topology().fields();
                std::string name = fields.id2name(s[0]);

                // switch the field in scene k, plot j
                window_.scene(k).primitive(j).set_active( name , rank );
                window_.scene(k).primitive(j).get_field_limits(ulim);
              }
              application_.write();
            }
            ImGui::SameLine();
            std::string cbar_label = "colorbar";
            static bool draw_colorbar = false;
            ImGui::Checkbox(cbar_label.c_str(), &draw_colorbar );
            if (draw_colorbar)
            {
              window_.draw_colorbar(application_.colormap(),ulim);
            }

            for (index_t i=0;i<entities.size();i++)
            {
              if (ImGui::TreeNode(entities[i].c_str()))
              {
                std::vector<std::string> entity = plot[entities[i]];
                for (index_t m=0;m<entity.size();m++)
                {
                  primitives_.insert( {entity[m],{k,j} } );

                  unsigned long address = std::stoul(entity[m],0,16);
                  Primitive* primitive = (Primitive*) address;

                  ImGui::Text("%s %lu",entities[i].c_str(),m);
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

                  if (primitive->number()==0) primitive->points_on() = primitive->triangles_on();
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
      ImGui::Combo("Mesh    ",&mesh_current, items.data() , items.size() );

      ImGui::SameLine();
      ImGui::SetNextItemWidth(5*width);
      static char iter_str[5] = "1";
      ImGui::InputText("Iter. ", iter_str , IM_ARRAYSIZE(iter_str));
      int nb_adapt = atoi(iter_str);

      ImGui::SetNextItemWidth(100);
      static int geometry_current = 0;
      const std::vector<std::string>& geometries = lib->geometries();
      items.resize( geometries.size() );
      for (index_t k=0;k<geometries.size();k++)
        items[k] = geometries[k].c_str();
      ImGui::Combo("Geometry",&geometry_current, items.data() , items.size() );

      ImGui::SameLine();
      static char output_name[128] = "output";
      ImGui::SetNextItemWidth(100);
      ImGui::InputText("Prefix", output_name, IM_ARRAYSIZE(output_name));

      ImGui::SetNextItemWidth(100);
      static int metric_current = 0;
      const std::vector<std::string>& metrics = lib->metrics();
      items.resize( metrics.size() );
      for (index_t k=0;k<metrics.size();k++)
        items[k] = metrics[k].c_str();
      ImGui::Combo("Metric  ",&metric_current, items.data() , items.size() );

      ImGui::SameLine();
      std::string mesh_name = meshes[mesh_current];
      if (ImGui::Button("Compute"))
      {
        ImGui::OpenPopup("adapt-mesh");
      }
      if (ImGui::BeginPopupModal("adapt-mesh", NULL, ImGuiWindowFlags_MenuBar))
      {
        ImGui::Text("Adapt mesh %s, metric %s, geometry %s?\n",mesh_name.c_str(),metrics[metric_current].c_str(),geometries[geometry_current].c_str());
        if (ImGui::Button("Adapt"))
        {

          ImGui::CloseCurrentPopup();
          std::string iter_cmd = "nb_iter="+std::to_string(nb_adapt);
          const char* command[] = {mesh_name.c_str(),geometries[geometry_current].c_str(),metrics[metric_current].c_str(),output_name,iter_cmd.c_str()};
          int result = programs::adapt(5,command);
          if (result!=0) printf("error in adaptation!\n");

          // check if this mesh needs to be updated
          printf("mesh name = %s\n",mesh_name.c_str());
          std::map<std::string,std::pair<index_t,index_t>>::iterator it;
          for (it=primitives_.begin();it!=primitives_.end();++it)
            printf("primitive %s -> scene = %lu, prim = %lu\n",it->first.c_str(),it->second.first,it->second.second);
          if (primitives_.find(mesh_name)!=primitives_.end())
          {
            // reload the mesh
            std::pair<index_t,index_t> plot = primitives_[mesh_name];
            index_t s = plot.first;  // scene
            index_t p = plot.second; // primitive

            // delete this root from the scene graph
            application_.remove( s , p );

            // restart the plotter
            application_.restart() = true;
          }
        }
        ImGui::SameLine();
        if (ImGui::Button("Cancel"))
        {
          ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
      }
    }

    if (ImGui::CollapsingHeader("Voronoi"))
    {
      std::vector<const char*> items;

      ImGui::SetNextItemWidth(100);
      static int mesh_current = 0;
      const std::vector<std::string>& meshes = lib->meshes();
      items.resize( meshes.size() );
      for (index_t k=0;k<meshes.size();k++)
        items[k] = meshes[k].c_str();
      ImGui::Combo("Mesh    ",&mesh_current, items.data() , items.size() );

      ImGui::SameLine();
      ImGui::SetNextItemWidth(5*width);
      static char iter_str[5] = "0";
      ImGui::InputText("Iter. ", iter_str , IM_ARRAYSIZE(iter_str));
      int nb_iter = atoi(iter_str);

      ImGui::SetNextItemWidth(100);
      static int geometry_current = 0;
      const std::vector<std::string>& geometries = lib->geometries();
      items.resize( geometries.size() + 1 );
      items[0] = "none";
      for (index_t k=0;k<geometries.size();k++)
        items[k+1] = geometries[k].c_str();
      ImGui::Combo("Geometry",&geometry_current, items.data() , items.size() );

      ImGui::SameLine();
      static bool hierarchical = false;
      ImGui::Checkbox("Hierarchical", &hierarchical );

      ImGui::SetNextItemWidth(100);
      static int sites_current = 0;
      const char* sites[] = {"points","sample","random","exact"};
      ImGui::Combo("Sites   ",&sites_current, sites , 3 );

      ImGui::SameLine();
      if (ImGui::Button("Calculate"))
      {
        // TODO treat CKF
        std::string mesh_name = meshes[mesh_current];

        printf("running voronoi with mesh %s, sites %s, geometry %s\n",mesh_name.c_str(),sites[sites_current],geometries[geometry_current+1].c_str());

        std::string iter_cmd = "nb_iter="+std::to_string(nb_iter);
        const char* command[] = {mesh_name.c_str(),sites[sites_current],geometries[geometry_current+1].c_str(),iter_cmd.c_str()};
        int result = programs::voronoi(4,command);
        UNUSED(result);
      }
    }

    if (ImGui::CollapsingHeader("Tools"))
    {
      ImGui::Text("Clipping:\n");

      ImGui::SetNextItemWidth(100);

      static float time = 0;
      ImGui::SliderFloat("time",&time,0,1,"%1.3f");

      std::string label = "modify";
      static bool mod = false;
      ImGui::Checkbox(label.c_str(), &mod );
      if (mod!=window_.modify_clipping_plane())
        window_.modify_clipping_plane() = mod;

      ImGui::SameLine();
      label = "show";
      ImGui::Checkbox(label.c_str(), &window_.show_clipping_plane() );

      ImGui::SameLine();
      if (ImGui::Button("flip normal"))
      {
        window_.flip_clipping_normal();
      }

      ImGui::SameLine();
      if (ImGui::Button("Clip"))
      {
        window_.clip();
      }
      if (ImGui::Button("Reset"))
      {
        window_.reset_clip();
      }

      ImGui::Separator();
      ImGui::Text("Axes:\n");
      label = "show axes";
      ImGui::Checkbox(label.c_str(), &window_.show_axes() );
      ImGui::SameLine();
      label = "center axes";
      ImGui::Checkbox(label.c_str(), &window_.center_axes() );

      ImGui::Separator();
      ImGui::Text("Colormap:");
      const char* colormaps[] = {"giraffe","parula","viridis","bgr","bwr","hsv","jet","hot"};
      static int current_colormap;
      ImGui::SetNextItemWidth(100);
      ImGui::Combo("Colormaps", &current_colormap, colormaps , 8 );
      ImGui::SameLine();
      if (ImGui::Button("change"))
      {
        application_.colormap().change_style(colormaps[current_colormap]);
        application_.write();
      }

      ImGui::Separator();
      ImGui::Text("Rendering:\n");
      static int fps = 60;
      ImGui::SetNextItemWidth(100);
      ImGui::SliderInt("fps",&fps,5,120,"%3d");
      window_.set_fps( fps );
    }

    if (ImGui::CollapsingHeader("Help"))
    {
      ImGui::Text("avro (c) Philip Caplan 2017-2020\nMiddlebury College, pcaplan@middlebury.edu\n");
      ImGui::Text("\n!!! warning !!!\nthis is a pre-alpha release so the interface is very rough!\n\n");
      ImGui::BulletText("rotate: hold ctrl (cmd) and click/move mouse");
      ImGui::BulletText("zoom: hold shift key and click/move mouse or scroll");
      ImGui::BulletText("pan: click/move mouse");
    }
  }

  ImGui::End();
}

void
Toolbar::end_draw() const
{
  ImGui::Render();
  #ifndef AVRO_HEADLESS_GRAPHICS
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
  #endif
}

#endif // AVRO_WITH_GL

} // graphics

} // avro
