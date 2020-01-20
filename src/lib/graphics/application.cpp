#include "graphics/application.h"

namespace avro
{

namespace graphics
{

void
ApplicationBase::write()
{
  // write all the data present in the scene graph
  printf("nb_scenes = %lu\n",scenes_.size());
  for (index_t k=0;k<scenes_.size();k++)
    scenes_[k]->write(manager_);
}

void
ShaderLibrary::set_matrices( SceneGraph& scene )
{
  // go through all the active shaders and assign the MVP and normalMatrix
  std::map<std::string,ShaderProgram>::iterator it;
  for (it=shaders_.begin();it!=shaders_.end();it++)
  {
    ShaderProgram& shader = it->second;
    shader.setUniform("MVP" , scene.mvp_matrix() );
    shader.setUniform("u_normalMatrix" , scene.normal_matrix() );
  }
}

void
OpenGL_Manager::write( Primitive& primitive )
{
  // allocate the buffers
  std::vector<GLuint> vbo(7);
  GL_CALL( glGenBuffers( vbo.size() , vbo.data() ) );
  GLuint& vbo_position  = vbo[0];
  GLuint& vbo_colour    = vbo[1];
  GLuint& vbo_normal    = vbo[2];
  GLuint& vbo_triangles = vbo[3];
  GLuint& vbo_edges     = vbo[4];
  GLuint& vbo_points    = vbo[5]; UNUSED(vbo_points);
  GLuint& vbo_feedback  = vbo[6];

  index_t nb_points = primitive.nb_points();

  if (primitive.number()>=2)
  {
    // bind the triangles
    std::vector<GLuint> triangles( primitive.triangles().begin() , primitive.triangles().end() );
    GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_triangles) );
    GL_CALL( glBufferData(GL_ELEMENT_ARRAY_BUFFER, triangles.size() * sizeof(GLuint), triangles.data() , GL_STATIC_DRAW) );
  }

  // bind the position buffer
  std::vector<GLfloat> points( primitive.points().begin() , primitive.points().end() );
  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, vbo_position) );
  GL_CALL( glBufferData(GL_ARRAY_BUFFER, points.size() * sizeof(GLfloat), points.data() , GL_STATIC_DRAW) );

  // bind the colour buffer
  std::vector<GLfloat> colors( primitive.colors().begin() , primitive.colors().end() );
  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, vbo_colour) );
  GL_CALL( glBufferData(GL_ARRAY_BUFFER, colors.size() * sizeof(GLfloat), colors.data() , GL_STATIC_DRAW) );

  // bind the normal buffer
  if (primitive.number()>=2)
  {
    std::vector<GLfloat> normals( primitive.normals().begin() , primitive.normals().end() );
    GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, vbo_normal) );
    GL_CALL( glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(GLfloat), normals.data() , GL_STATIC_DRAW) );
  }

  // bind the triangle data
  GLuint id_vao_triangles;
  GL_CALL( glGenVertexArrays( 1, &id_vao_triangles ) );
  GL_CALL( glBindVertexArray(id_vao_triangles) );
  vao_triangles_.insert( {&primitive,id_vao_triangles} );

  GL_CALL( glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, vbo_triangles ) );

  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, vbo_position ) );
  GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
  GL_CALL( glEnableVertexAttribArray(0) );

  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, vbo_colour ) );
  GL_CALL( glVertexAttribPointer( 1, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
  GL_CALL( glEnableVertexAttribArray(1) );

  if (primitive.number()>=2)
  {
    GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, vbo_normal ) );
    GL_CALL( glVertexAttribPointer( 2, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
    GL_CALL( glEnableVertexAttribArray(2) );
  }

  // bind the edges
  GLuint id_vao_edges;
  GL_CALL( glGenVertexArrays( 1, &id_vao_edges ) );
  GL_CALL( glBindVertexArray(id_vao_edges) );
  vao_edges_.insert( {&primitive,id_vao_edges} );

  std::vector<GLuint> edges( primitive.edges().begin() , primitive.edges().end() );
  GL_CALL( glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, vbo_edges ) );
  GL_CALL( glBufferData(GL_ELEMENT_ARRAY_BUFFER, edges.size() * sizeof(GLuint), edges.data() , GL_STATIC_DRAW) );

  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, vbo_position ) );
  GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
  GL_CALL( glEnableVertexAttribArray(0) );

  // bind the points
  GLuint id_vao_points;
  GL_CALL( glGenVertexArrays( 1, &id_vao_points ) );
  GL_CALL( glBindVertexArray(id_vao_points) );
  vao_points_.insert( {&primitive,id_vao_points} );

  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, vbo_position ) );
  GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
  GL_CALL( glEnableVertexAttribArray(0) );

  // reset the vao bound to the gl
  GL_CALL( glBindVertexArray(0) );
}

void
OpenGL_Manager::draw( Primitive& primitive , TransformFeedbackResult* feedback )
{
  // indicate to the gl that we want to use the shader
  ShaderProgram shader = *shader_.at(&primitive);
  shader.use();

  if (feedback!=nullptr)
  {
    GLuint& query = feedback->query();
    GLuint& buffer = feedback->buffer();

    glGenQueries(1, &query );
    GL_CALL( glBeginQuery( GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN , query ) );

    GL_CALL( glGenBuffers( 1 , &buffer ) );
    GL_CALL( glBindBuffer( GL_TRANSFORM_FEEDBACK_BUFFER , buffer ) );
    GL_CALL( glBufferData( GL_TRANSFORM_FEEDBACK_BUFFER , 8*primitive.nb_triangles()*sizeof(GLfloat) , NULL , GL_DYNAMIC_COPY ) );
    GL_CALL( glBindBufferBase( GL_TRANSFORM_FEEDBACK_BUFFER , 0 , buffer ) );
    GL_CALL( glBindVertexArray( feedback->vao() ) );

    GL_CALL( glEnable(GL_RASTERIZER_DISCARD) );
    GL_CALL( glBeginTransformFeedback(GL_TRIANGLES) );
  }

  // bind the vao associated with this primitive
  if (primitive.number()>=2)
  {
    // draw the triangles
    if (primitive.triangles_on())
    {
      GL_CALL( glBindVertexArray( vao_triangles_.at(&primitive) ) );
      GL_CALL( glDrawElements(GL_TRIANGLES, 3*primitive.nb_triangles() , GL_UNSIGNED_INT , 0 ) )
    }
  }

  if (feedback!=nullptr)
  {
    GL_CALL( glEndTransformFeedback() );
    GL_CALL( glDisable(GL_RASTERIZER_DISCARD) );

    GLuint nb_primitives = 0;
    GL_CALL( glEndQuery( GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN ) );
    GL_CALL( glGetQueryObjectuiv( feedback->query() , GL_QUERY_RESULT , &nb_primitives ) );
    printf("primitives written = %u\n",nb_primitives);

    // there are 4 coordinates for every output and 3 vertices per prim (color + coord)
    index_t nb_data_per_prim = 4*3*2;
    index_t size = nb_primitives*nb_data_per_prim;

    GLfloat* feedbackBuffer = (GLfloat*) malloc( size*sizeof(GLfloat) );
    GL_CALL( glBindBuffer( GL_TRANSFORM_FEEDBACK_BUFFER , feedback->buffer() ) );
    GL_CALL( glGetBufferSubData( GL_TRANSFORM_FEEDBACK_BUFFER , 0 , size*sizeof(GLfloat) , feedbackBuffer ) );

    //printBuffer(feedbackBuffer, nb_primitives , nb_data_per_prim );
  }

  if (primitive.number()>=1)
  {
    if (primitive.edges_on())
    {
      GL_CALL( glBindVertexArray( vao_edges_.at(&primitive) ) );
      GL_CALL( glDrawElements( GL_LINES , primitive.nb_edges()*2 , GL_UNSIGNED_INT , 0 ) )
    }
  }

  // draw the points
  if (primitive.points_on())
  {
    GL_CALL( glBindVertexArray(vao_points_.at(&primitive) ) );
    GL_CALL( glPointSize(10.0f) );
    GL_CALL( glDrawArrays( GL_POINTS , 0 , primitive.nb_points() ) );
  }

  // reset the vao bound to the gl
  GL_CALL( glBindVertexArray(0) );
}

void
OpenGL_Manager::draw( SceneGraph& scene , TransformFeedbackResult* feedback )
{
  if (!scene.update()) return;

  shaders_.set_matrices(scene);

  for (index_t k=0;k<scene.nb_primitives();k++)
  {
    draw( scene.primitive(k) , feedback );
  }

  scene.set_update(false);
}

static Trackball* trackball_instance = nullptr;

static void
mouse_button_callback(GLFWwindow* win ,int button,int action,int mods)
{
  if (action == GLFW_PRESS)
  {
    double xpos,ypos;
    glfwGetCursorPos(win,&xpos,&ypos);
    trackball_instance->MouseDown(button,action,mods,(int)xpos,(int)ypos);
  }
  if (action == GLFW_RELEASE)
  {
    trackball_instance->MouseUp();
  }
}

static void
mouse_move_callback(GLFWwindow* win, double xpos, double ypos)
{
  trackball_instance->MouseMove((int)xpos,(int)ypos);
}

static void
mouse_scroll_callback(GLFWwindow* win, double xpos, double ypos)
{
  trackball_instance->MouseWheel(xpos,ypos);
}

static void
keyboard_callback(GLFWwindow* win,int key,int scancode,int action,int mods)
{
  if (action == GLFW_PRESS)
  {
    trackball_instance->KeyDown(key);
  }

  if (action == GLFW_RELEASE)
  {
    trackball_instance->KeyUp();
  }
}

GLFW_Window::GLFW_Window( int width , int height ) :
  title_("avro window"),
  width_(width),
  height_(height),
  camera_(glm::vec3(0.0f,0.0f,7.0f)),
  trackball_(&camera_,glm::vec4(0.0f,0.0f,(float)width_,(float)height_))
{
  window_ = glfwCreateWindow( width_ , height_ , title_.c_str() , NULL, NULL);
  if (!window_)
  {
    glfwTerminate();
    avro_assert_not_reached;
  }
  glfwMakeContextCurrent(window_);

  trackball_instance = &trackball_;

  // initialize the trackball callbacks
  glfwSetCursorPosCallback(window_,&mouse_move_callback);
  glfwSetMouseButtonCallback(window_,&mouse_button_callback);
  glfwSetScrollCallback(window_,&mouse_scroll_callback);
  glfwSetKeyCallback(window_,&keyboard_callback);
}

template class Application<GLFW_Interface<OpenGL_Manager>>;
template class Application<GLFW_Interface<Vulkan_Manager>>;

} // graphics

} // avro
