#include "graphics/gui.h"
#include "graphics/primitives.h"
#include "graphics/vao.h"
#include "graphics/window.h"

#include "graphics/gl.h"

#include <imgui/GL/imgui_impl_glfw.h>
#include <imgui/GL/imgui_impl_opengl3.h>

#define EMPTY_LABEL empty_label().c_str()
#define UNIQUE_LABEL(X) unique_label(X).c_str()

namespace avro
{

namespace graphics
{

index_t label_counter;

GUI::GUI( Window& window ) :
 window_(window),
 count_(5)
{
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

std::string
empty_label() {
  return "##" + std::to_string(label_counter++);
}

std::string
unique_label( const std::string& s ) {
  return s + "##" + std::to_string(label_counter++);
}

static std::vector<const char*>
convert_to_char( const std::vector<std::string>& s ) {
  std::vector<const char*> cs( s.size() );
  for (index_t k = 0; k < s.size(); k++)
    cs[k] = s[k].c_str();
  return cs;
}

void
GUI::draw() {

  //bool capture_mouse = ImGui::GetIO().WantCaptureMouse;
  ImGuiWindowFlags window_flags = 0;
  window_flags |= ImGuiWindowFlags_NoCollapse;
  window_flags |= ImGuiWindowFlags_NoResize;
  //window_flags |= ImGuiWindowFlags_AlwaysHorizontalScrollbar;

  bool capture_mouse = ImGui::GetIO().WantCaptureMouse;
  window_.needs_drawing(true);
  window_.enable_controls( !capture_mouse );

  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();

  // set the controls at the top-left of the window
  ImGui::SetNextWindowPos( ImVec2( 0 , 0 ) );
  ImGui::SetNextWindowSize( ImVec2( 256 , window_.height() ) );

  std::vector<std::string> entities = {"Faces","Edges"};//,"Nodes"};

  bool active = true;
  ImGui::SetNextItemWidth(200);

  label_counter = 0;

  if (ImGui::Begin("Controls",&active,window_flags))
  {
    if (ImGui::Button("PNG")) {
      ImGui::Text("enter file name");
      ImGui::SameLine();
      ImGui::Button("save file!");

    }
    ImGui::SameLine();
    if (ImGui::Button("EPS")) {

      ImGui::OpenPopup("eps-save");
    }

    if (ImGui::BeginPopup("eps-save")) {

      ImGui::Button("save");
      ImGui::EndPopup();
    }

    ImGui::SameLine();
    if (ImGui::Button("View")) {
      printf("load view matrix (or parameters)\n");
    }

    ImGui::SameLine();
    ImGui::Checkbox("lighting",&window_.lighting());

    ImGui::Separator();

    ImGui::Text("Plots");
    ImGui::SameLine();
    if (ImGui::Button("Hide All")) {
      for (index_t i = 0; i < window_.nb_plots(); i++) {
        window_.plot(i).hidden() = true;
        for (index_t k = 0; k < window_.plot(i).active_vao().nb_triangles(); k++)
          window_.plot(i).active_vao().triangles(k).visible() = false;
        for (index_t k = 0; k < window_.plot(i).active_vao().nb_edges(); k++)
          window_.plot(i).active_vao().edges(k).visible() = false;
      }
    }

    ImGui::SameLine();
    ImGui::SetNextItemWidth(100);
    static int colormap_index = 0;
    const char* colormaps[] = {"[colormap]","viridis","parula","giraffe","bwr","hsv","jet"};
    if (ImGui::Combo(empty_label().c_str(),&colormap_index,colormaps,7)) {
      if (colormap_index == 0) colormap_index = 1;
      window_.select_colormap(colormaps[colormap_index]);
    }
    ImGui::Separator();

    std::string label;
    for (index_t i = 0; i < window_.nb_plots(); i++) {

      // get the active vao information
      const nlohmann::json& info = window_.plot(i).active_vao().get_info();

      std::string label = "Plot" + std::to_string(i);
      if (ImGui::TreeNode(label.c_str())) {

        static int current_mesh = 1;

        std::vector<std::string> vao_labels = window_.plot(i).vao_labels();
        vao_labels.insert( vao_labels.begin() , "[group]");
        std::vector<const char*> c_vao_labels = convert_to_char(vao_labels);

        ImGui::SetNextItemWidth(100);
        label = empty_label();
        if (ImGui::Combo(label.c_str(),&current_mesh,c_vao_labels.data(),c_vao_labels.size())) {
          // pick which vao we need to render
          if (current_mesh == 0) current_mesh = 1;
          window_.plot(i).set_active( current_mesh-1 );
        }

        VertexAttributeObject& vao = window_.plot(i).active_vao();

        ImGui::SameLine();
        if (!window_.plot(i).hidden()) {
          if (ImGui::Button("Hide")) {
            // completely turn off the plot
            for (index_t k = 0; k < vao.nb_triangles(); k++)
              vao.triangles(k).visible() = false;
            for (index_t k = 0; k < vao.nb_edges(); k++)
              vao.edges(k).visible() = false;
            window_.plot(i).hidden() = true;
          }
        }
        else {
          if (ImGui::Button("Show")) {
            // completely turn on the plot
            for (index_t k = 0; k < vao.nb_triangles(); k++)
              vao.triangles(k).visible() = true;
            for (index_t k = 0; k < vao.nb_edges(); k++)
              vao.edges(k).visible() = true;
            window_.plot(i).hidden() = false;
          }
        }

        int current_field = 1;
        int current_field_rank = 1;
        std::vector<nlohmann::json> fields = info["fields"];

        std::vector<std::string> field_names = info["field_names"];
        std::vector<std::string> field_ranks;
        std::vector<int> rank_index;
        if (field_names.size() > 0) {
          nlohmann::json jf = fields[0];
          std::vector<std::string> tmps = jf.at("ranks");
          field_ranks.assign( tmps.begin() , tmps.end() );

          std::vector<int> tmpi = jf.at("rank_index");
          rank_index.assign(tmpi.begin(),tmpi.end());
        }

        field_names.insert( field_names.begin() , "geometry" );
        field_names.insert( field_names.begin() , "constant" );
        field_names.insert( field_names.begin() , "[field]" );
        field_ranks.insert( field_ranks.begin() , "[rank]" );

        std::vector<const char*> field_labels = convert_to_char(field_names);
        std::vector<const char*> rank_labels = convert_to_char(field_ranks);

        ImGui::SetNextItemWidth(100);
        label = empty_label();
        if (ImGui::Combo(label.c_str(),&current_field,field_labels.data(),field_labels.size())) {

          // prevent from picking the label
          if (current_field == 0) current_field = 1;

          if (current_field == 1) {
            vao.geometry_color() = false;
            vao.uniform_color()  = true;
          }
          else if (current_field == 2) {
            vao.geometry_color() = true;
            vao.uniform_color()  = false;
          }
          else if (vao.nb_fields() > 0) {
            vao.geometry_color() = false;
            vao.uniform_color()  = false;
            vao.set_field(field_names[current_field]);

            nlohmann::json jf = fields[current_field-3];
            std::vector<std::string> tmps = jf.at("ranks");
            field_ranks.assign( tmps.begin() , tmps.end() );

            std::vector<int> tmpi = jf.at("rank_index");
            rank_index.assign(tmpi.begin(),tmpi.end());
          }
        }
        ImGui::SameLine();
        ImGui::SetNextItemWidth(80);
        label = empty_label();
        if (ImGui::Combo(label.c_str(),&current_field_rank,rank_labels.data(),rank_labels.size())) {
          if (vao.nb_fields() > 0) {
            if (current_field_rank == 0) current_field_rank = 1;
            if (!vao.geometry_color() && !vao.uniform_color())
              vao.set_rank(rank_index[current_field_rank-1]);
          }
        }

        label = unique_label("min");
        ImGui::SetNextItemWidth(75);
        ImGui::InputFloat(label.c_str(),&window_.plot(i).active_vao().umin(),-1e20,1e20,"%.1f");

        ImGui::SameLine();
        label = unique_label("max");
        ImGui::SetNextItemWidth(75);
        ImGui::InputFloat(label.c_str(),&window_.plot(i).active_vao().umax(),-1e20,1e20,"%.1f");

        ImGui::Separator();
        ImGui::SetNextItemWidth(50);
        label = unique_label("tessellation");
        ImGui::SliderInt(label.c_str(),&vao.tessellation_level(),1,20);

        ImGui::SetNextItemWidth(80);
        const char* clip_styles[] = {"[off]","pixel","prim."};
        label = unique_label("clipping");
        if (ImGui::Combo(label.c_str(),&window_.plot(i).clip().style(),clip_styles,3)) {
        }
        ImGui::SameLine();
        label = unique_label("show");
        if (ImGui::Checkbox(label.c_str(),&window_.plot(i).clip().visible())) {
          window_.plot(i).clip().update();
        }
        const char* clip_dims[] = {"X","Y","Z"};
        label = empty_label();
        ImGui::SetNextItemWidth(80);
        if (ImGui::Combo(label.c_str(),&window_.plot(i).clip().dimension(),clip_dims,3)) {}
        ImGui::SameLine();
        label = unique_label("flip");
        if (ImGui::Button(label.c_str())) {
          window_.plot(i).clip().flip();
        }

        // slider to control the clip plane distance
        ImGui::SetNextItemWidth(100);
        label = unique_label("distance");
        float length_scale = window_.plot(i).clip().length_scale();
        ImGui::SliderFloat(label.c_str(),&window_.plot(i).clip().distance(),-length_scale,length_scale);
        ImGui::Separator();

        // update the clip
        window_.plot(i).clip().update();


        if (ImGui::TreeNode("Faces")) {

          // get all the entities from the plot vao
          std::vector<nlohmann::json> facets = info["triangles"];
          for (index_t k = 0; k < facets.size(); k++) {

            std::string s = facets[k]["name"];
            label = unique_label(s);
            if (ImGui::Checkbox(label.c_str(),&vao.triangles(k).visible())) {
            }
          }
          ImGui::TreePop();
        }

        if (ImGui::TreeNode("Edges")) {

          // get all the entities from the plot vao
          std::vector<nlohmann::json> facets = info["edges"];
          for (index_t k = 0; k < facets.size(); k++) {

            std::string s = facets[k]["name"];
            label = unique_label(s);
            if (ImGui::Checkbox(label.c_str(),&vao.edges(k).visible())) {
            }
          }
          ImGui::TreePop();
        }

        ImGui::TreePop();
      }
    }
    ImGui::Separator();

    real_t memory = 0;
    for (index_t k = 0; k < window_.nb_plots(); k++) {
      for (index_t j = 0; j < window_.plot(k).nb_vao(); j++)
        memory += window_.plot(k).vao(j).get_memory();
    }

    ImGui::Text("FPS %.1f", ImGui::GetIO().Framerate);
    ImGui::Text("Draw count: %lu" , window_.draw_count());

    if (memory < 1e6)
      ImGui::Text("GPU mem. in use: %.1f kB" ,real_t(memory)/1e3);
    else
      ImGui::Text("GPU mem. in use: %.1f MB" ,real_t(memory)/1e6);

    ImGui::Separator();

    ImGui::End();
  }

  ImGui::Render();
  window_.draw(false); // do not swap buffers, we need to render the gui first
  #if AVRO_HEADLESS_GRAPHICS == 0
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
  #endif

  glfwSwapBuffers(window_.window());
}


} // graphics

} // avro
