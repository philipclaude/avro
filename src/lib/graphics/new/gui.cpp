#include "graphics/new/gui.h"
#include "graphics/new/window.h"

#include "graphics/gl.h"

#include <imgui/GL/imgui_impl_glfw.h>
#include <imgui/GL/imgui_impl_opengl3.h>

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
header( const std::string& name) {
  return "[" + name + "]";
}

static std::vector<const char*>
convert_to_char( const std::vector<std::string>& s , const std::string& H = std::string() ) {
  std::vector<const char*> cs( s.size() );

  index_t idx = 0;
  if (!H.empty()) {
    cs.insert(cs.begin() , H.c_str() );
    idx = 1;
  }

  for (index_t k = 0; k < s.size(); k++)
    cs[idx++] = s[k].c_str();
  return cs;
}

void
GUI::draw() {

  //bool capture_mouse = ImGui::GetIO().WantCaptureMouse;
  ImGuiWindowFlags window_flags = 0;
  //window_flags |= ImGuiWindowFlags_NoCollapse;
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
  ImGui::SetNextWindowSize( ImVec2( 200 , window_.height() ) );

  std::vector<std::string> entities = {"Faces","Edges","Nodes"};

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
    ImGui::Separator();

    ImGui::Text("Plots");
    ImGui::SameLine();
    if (ImGui::Button("Hide All")) {
      printf("hide all plots\n");
    }
    ImGui::Separator();

    for (index_t i = 0; i < window_.nb_plots(); i++) {

      // get the active vao information
      const nlohmann::json& info = window_.plot(i).active_vao().get_info();
      //std::cout << info.dump() << std::endl;

      std::string label = "Plot" + std::to_string(i);
      if (ImGui::TreeNode(label.c_str())) {

        static int current_mesh = 1;

        const std::vector<std::string>& vao_labels = window_.plot(i).vao_labels();
        std::vector<const char*> c_vao_labels = convert_to_char(vao_labels,header("mesh"));

        for (index_t k = 0; k < c_vao_labels.size(); k++)
          printf("label %lu = %s\n",k,c_vao_labels[k]);

        ImGui::SetNextItemWidth(100);
        if (ImGui::Combo(empty_label().c_str(),&current_mesh,c_vao_labels.data(),c_vao_labels.size())) {
          // pick which vao we need to render
        }
        ImGui::SameLine();
        if (ImGui::Button("Hide")) {
          // completely turn off the plot
        }

        static int current_field = 0;
        static int current_field_rank = 0;
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

        std::vector<const char*> field_labels = convert_to_char(field_names,header("field"));
        std::vector<const char*> rank_labels = convert_to_char(field_ranks,header("rank"));

        ImGui::SetNextItemWidth(100);
        if (ImGui::Combo("##1",&current_field,field_labels.data(),field_labels.size())) {

        }
        ImGui::SameLine();
        ImGui::SetNextItemWidth(100);
        if (ImGui::Combo("##2",&current_field_rank,rank_labels.data(),rank_labels.size())) {

        }


        ImGui::Separator();
        static int level = 1;
        ImGui::SetNextItemWidth(50);
        ImGui::SliderInt("tessellation##1",&level,1,20);

        static bool show_clip = false;
        static bool modify_clip = false;

        ImGui::SetNextItemWidth(100);
        const char* clip_styles[] = {"[off]","pixel","primitive"};
        static int current_clip_style = 0;
        if (ImGui::Combo("clipping##3",&current_clip_style,clip_styles,3)) {

        }
        if (ImGui::Checkbox("show##1",&show_clip)) {

        }
        ImGui::SameLine();
        if (ImGui::Checkbox("modify##1",&modify_clip)) {

        }
        ImGui::Separator();

        for (index_t j = 0; j < entities.size(); j++) {


          if (ImGui::TreeNode(entities[j].c_str())) {

            // get all the entities from the plot vao

            ImGui::TreePop();
          }
        }

        ImGui::TreePop();
      }
    }
    ImGui::Separator();

    ImGui::Text("FPS %.1f", ImGui::GetIO().Framerate);
    ImGui::Text("Draw count: %lu" , window_.draw_count());

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
