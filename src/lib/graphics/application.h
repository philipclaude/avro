#ifndef avro_LIB_GRAPHICS_CLIENT_H_
#define avro_LIB_GRAPHICS_CLIENT_H_

#include "common/error.h"
#include "common/tools.h"
#include "common/types.h"

#include "graphics/controls.h"
#include "graphics/primitive.h"
#include "graphics/shader.h"

#include "library/eps.h"

#include <map>
#include <memory>
#include <unordered_set>
#include <string>
#include <vector>

namespace avro
{

class Mesh;
class TopologyBase;

namespace graphics
{

class SceneGraph
{
public:
  typedef std::shared_ptr<Primitive> Primitive_ptr;

  SceneGraph() :
    update_(true)
  {}

  index_t add_primitive( TopologyBase& topology )
  {
    index_t id = primitive_.size();
    Primitive_ptr primitive = std::make_shared<Primitive>(topology,this);
    primitive_.push_back(primitive);
    return id;
  }

  void write( GraphicsManager& manager )
  {
    for (index_t k=0;k<primitive_.size();k++)
      primitive_[k]->write(manager);
  }

  bool update() const { return update_; }
  void set_update( bool x ) { update_ = x; }

  void update_matrices( const Trackball& trackball , float,float,float );

  const mat4& mvp_matrix() const { return mvp_matrix_; }
  const mat4& normal_matrix() const { return normal_matrix_; }

  mat4& mvp_matrix() { return mvp_matrix_; }
  mat4& normal_matrix() { return normal_matrix_; }

  index_t nb_primitives() const { return primitive_.size(); }

  Primitive& primitive( index_t k ) { return *primitive_[k].get(); }

private:
  std::vector<Primitive_ptr> primitive_; // roots of the scene graph

  // store all the matrices here
  mat4 mvp_matrix_;
  mat4 view_matrix_;
  mat4 proj_matrix_;
  mat4 model_matrix_;
  mat4 normal_matrix_;

  bool update_;
};

class TransformFeedbackResult
{
public:

  void append_triangles( const GLfloat* buffer , index_t nb_triangles )
  {
    index_t n = 0;
    for (index_t k=0;k<nb_triangles;k++)
    {
      //printf("primitve %lu\n",k);
      for (index_t j=0;j<3;j++)
      {
        for (coord_t d=0;d<3;d++)
          triangle_points_.push_back( buffer[n+d] );
        for (coord_t d=0;d<3;d++)
          triangle_colors_.push_back( buffer[n+4+d] );

        /*printf("v = (%g,%g,%g,%g), c = (%g,%g,%g,%g)\n",
                buffer[n  ],buffer[n+1],buffer[n+2],buffer[n+3] ,
                buffer[n+4],buffer[n+5],buffer[n+6],buffer[n+7]);*/
        n += 8;
      }
    }
  }

  void append_edges( const GLfloat* buffer , index_t nb_edges )
  {
    index_t n = 0;
    for (index_t k=0;k<nb_edges;k++)
    {
      //printf("primitve %lu\n",k);
      for (index_t j=0;j<2;j++)
      {
        for (coord_t d=0;d<3;d++)
          edge_points_.push_back( buffer[n+d] );
        for (coord_t d=0;d<3;d++)
          edge_colors_.push_back( buffer[n+4+d] );

        /*printf("v = (%g,%g,%g,%g), c = (%g,%g,%g,%g)\n",
                buffer[n  ],buffer[n+1],buffer[n+2],buffer[n+3] ,
                buffer[n+4],buffer[n+5],buffer[n+6],buffer[n+7]);*/
        n += 8;
      }
    }
  }

  const std::vector<real_t>& triangle_points() const { return triangle_points_; }
  const std::vector<real_t>& triangle_colors() const { return triangle_colors_; }

  const std::vector<real_t>& edge_points() const { return edge_points_; }
  const std::vector<real_t>& edge_colors() const { return edge_colors_; }

private:
  std::vector<real_t> triangle_points_;
  std::vector<real_t> triangle_colors_;

  std::vector<real_t> edge_points_;
  std::vector<real_t> edge_colors_;

  std::vector<real_t> vertex_points_;
  std::vector<real_t> vertex_colors_;
};

class GraphicsManager
{
public:
  virtual ~GraphicsManager() {}
  virtual void write( Primitive& primitive ) = 0;
  virtual void draw( SceneGraph& scene , TransformFeedbackResult* feedback=nullptr ) = 0;
};

class GLFW_Window
{
public:
  GLFW_Window( GraphicsManager& manager , int width , int height , const std::string& title );

  void
  mouse_button_callback(int button,int action,int mods)
  {
    if (action == GLFW_PRESS)
    {
      double xpos,ypos;
      glfwGetCursorPos(window_,&xpos,&ypos);
      trackball_.MouseDown(button,action,mods,(int)xpos,(int)ypos);
    }
    if (action == GLFW_RELEASE)
    {
      trackball_.MouseUp();
    }
  }

  void
  mouse_move_callback(double xpos, double ypos)
  {
    trackball_.MouseMove((int)xpos,(int)ypos);
  }

  void
  mouse_scroll_callback(double xpos, double ypos)
  {
    trackball_.MouseWheel(xpos,ypos);
  }

  void key_callback(int key, int scancode, int action, int mods)
  {
    trackball_.KeyDown(key);
  }

  void make_current()
  {
    glfwMakeContextCurrent(window_);
  }

  void update_view();

  bool should_close()
  {
    return glfwWindowShouldClose(window_) || (glfwGetKey(window_, GLFW_KEY_ESCAPE ) == GLFW_PRESS);
  }

  void save_eps( const std::string& filename )
  {
    TransformFeedbackResult feedback;
    for (index_t k=0;k<scene_.size();k++)
    {
      printf("capture with transform feedback!!\n");
      manager_.draw(scene_[k],&feedback);
    }


    library::epsFile eps;
    int viewport[4] = {0,0,1024,640};
    eps.set_viewport(viewport);
    eps.add_triangles( feedback.triangle_points() , feedback.triangle_colors() );
    eps.add_edges( feedback.edge_points() , feedback.edge_colors() );
    eps.write( filename );
  }

  void setup()
  {
    // ensure we can capture the escape key being pressed below
    glfwSetInputMode(window_, GLFW_STICKY_KEYS, GL_TRUE);

    // hide the mouse and enable unlimited mouvement
    //glfwSetInputMode(window_, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    glfwPollEvents();
    glfwSetCursorPos(window_, width_/2, height_/2);

    glClearColor(0.0f, 0.0f, 0.4f, 0.0f);

    // enable depth test
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);

    // option to cull triangles which normal is not towards the camera
    glDisable(GL_CULL_FACE);

    //glEnable(GL_BLEND);
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0, 1.5);
  }

  void begin_draw()
  {
    make_current();

    glfwPollEvents();

    glClearColor (1.0, 1.0, 1.0, 0.0); // white
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //update_view();
  }

  void end_draw()
  {
    glfwSwapBuffers(window_);
  }

  index_t nb_scene() const { return scene_.size(); }
  SceneGraph& scene( index_t k ) { return scene_[k]; }

  index_t create_scene()
  {
    index_t id = scene_.size();
    scene_.push_back(SceneGraph());
    return id;
  }

private:
  std::string title_;
  GLFWwindow* window_;
  GraphicsManager& manager_;

  int width_;
  int height_;
  float fov_ = 45.0f;

  vec3 position_;
  float angles_[2];

  std::vector<SceneGraph> scene_;
  Trackball trackball_;
  Camera camera_;
};

class ShaderLibrary
{
public:
  ShaderLibrary()
  {
    // generate all the shaders!!
  }

  const ShaderProgram& operator[] ( const std::string& name ) const
  {
    avro_assert( shaders_.find(name)!=shaders_.end() );
    return shaders_.at(name);
  }

  ShaderProgram& operator[] ( const std::string& name )
  {
    avro_assert( shaders_.find(name)!=shaders_.end() );
    return shaders_.at(name);
  }

  void set_matrices( SceneGraph& scene );

  void create();

private:

  // store all the shaders
  std::map<std::string,ShaderProgram> shaders_;
};

class OpenGL_Manager : public GraphicsManager
{
  template<typename type> friend class Application;

private:
  OpenGL_Manager()
  {}

public:
  void write( Primitive& primitive );

  void draw( GLFW_Window& window , TransformFeedbackResult* feedback=nullptr );
  void draw( SceneGraph& scene , TransformFeedbackResult* feedback=nullptr );

  void select_shader( Primitive& primitive , const std::string& name )
  { shader_[&primitive] = &shaders_[name]; }
  void create_shaders() { shaders_.create(); }

private:

  void set_matrices( SceneGraph& scene );
  void draw( Primitive& primitive , TransformFeedbackResult* feedback=nullptr );

  // map from primitive to vao
  std::map<Primitive*,index_t> vao_points_;
  std::map<Primitive*,index_t> vao_edges_;
  std::map<Primitive*,index_t> vao_triangles_;

  std::map<Primitive*,index_t> vao_feedback_points_;
  std::map<Primitive*,index_t> vao_feedback_edges_;
  std::map<Primitive*,index_t> vao_feedback_triangles_;

  std::map<Primitive*,ShaderProgram*> shader_;
  ShaderLibrary shaders_;
};

class Vulkan_Manager : public GraphicsManager
{
  template<typename type> friend class Application;

private:
  Vulkan_Manager() {}

  void write( Primitive& primitive ) { avro_implement; }
  void draw( SceneGraph& scene , TransformFeedbackResult* feedback=nullptr )
  { avro_implement; }

  void create_shaders()
  { avro_implement; }

};

class WV_Manager : public GraphicsManager
{
  template<typename type> friend class Application;

private:
  WV_Manager();

  void write( Primitive& primitive ) { avro_implement; }
  void draw( SceneGraph& scene , TransformFeedbackResult* feedback=nullptr )
  { avro_assert_not_reached; }

};

class ApplicationBase
{
protected:
  ApplicationBase( GraphicsManager& manager ) :
    manager_(manager)
  {}

  virtual ~ApplicationBase() {}

  virtual void run() = 0;
  void write();

protected:
  std::vector<SceneGraph*> scenes_;

private:
  GraphicsManager& manager_;
};

template<typename type> class Application;

template<typename T> struct GLFW_Interface;
struct Web_Interface;

template<typename API_t>
class Application<GLFW_Interface<API_t>> : public ApplicationBase
{
public:
  Application() :
    ApplicationBase(manager_)
  {
    // initialize OpenGL
    avro_assert_msg( glfwInit() , "problem initializing OpenGL!" );

    // set the version
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  }

  void initialize()
  {
    avro_assert( window_.size()!=0 );
    window_[0]->make_current();

    // load gl stuff and print info
    gladLoadGL();
    dumpGLInfo();

    manager_.create_shaders();
  }

  void run()
  {
    for (index_t k=0;k<window_.size();k++)
      window_[k]->setup();

    write();

    for (index_t k=0;k<window_.size();k++)
      window_[k]->update_view();

    const double fps = 1.0 / 60.0;
    double last_update_time = 0;  // number of seconds since the last loop
    double last_frame_time = 0;   // number of seconds since the last frame

     // start the rendering loop
     bool done = false;
     while (!done)
     {

       double now = glfwGetTime();

       // This if-statement only executes once every 60th of a second
       if ((now - last_frame_time) >= fps)
       {
         // draw your frame here
         for (index_t k=0;k<window_.size();k++)
         {
           window_[k]->begin_draw();

           for (index_t j=0;j<window_[k]->nb_scene();j++)
             manager_.draw(window_[k]->scene(j));

           window_[k]->end_draw();

           if (window_[k]->should_close())
           {
             printf("window %lu requested close\n",k);
             done = true;
           }

         }

         // only set lastFrameTime when you actually draw something
         last_frame_time = now;
        }

        // set lastUpdateTime every iteration
        last_update_time = now;
     }
  }

protected:
  void add_window( GLFW_Window* window )
    { window_.push_back(window); }

protected:
  API_t manager_;
  std::vector<GLFW_Window*> window_;
};

template<>
class Application<Web_Interface> : public ApplicationBase
{
public:
  Application() :
    ApplicationBase(manager_)
  {}

  void run() { avro_implement; } // run the server

  void save_eps() // will also need transformation matrices sent from the web
  {
    // first set the transformation to the scene
    OpenGL_Manager manager_gl;
    scene_.write(manager_gl);
    TransformFeedbackResult feedback;
    manager_gl.draw(scene_,&feedback);
  }

protected:
  WV_Manager manager_;
  SceneGraph scene_;
};

class Visualizer : public Application<GLFW_Interface<OpenGL_Manager>>
{
public:
  Visualizer() :
    main_(manager_,1024,612,"3D")
    //side_(manager_,512,306,"2D")
  {
    add_window( &main_ );
    //add_window( &side_ );

    initialize();
  }

  void add_topology( TopologyBase& topology )
  {
    // add the topology to the relevant windows
    index_t id = main_.create_scene();
    scenes_.push_back( &main_.scene(id) );
    index_t prim_id = main_.scene(id).add_primitive(topology); // create a new root in the scene graph
    manager_.select_shader( main_.scene(id).primitive(prim_id) , "wv" );

    #if 0
    // add the topology to the relevant windows
    id = side_.create_scene();
    scenes_.push_back( &side_.scene(id) );
    prim_id = side_.scene(id).add_primitive(topology); // create a new root in the scene graph
    manager_.select_shader( side_.scene(id).primitive(prim_id) , "wv" );
    #endif
  }

  GLFW_Window main_;
  //GLFW_Window side_;
};

class WebVisualizer : public Application<Web_Interface>
{
public:

  void add_topology( TopologyBase& topology )
  {
    // create a new root in the scene graph
    scene_.add_primitive(topology);
  }

private:
  SceneGraph scene_;

};

} // graphics

} // avro

#endif
