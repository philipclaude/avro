#include "unit_tester.hpp"

#include "graphics/gl.h"
#include "graphics/math.h"

using namespace avro;
using namespace avro::graphics;

UT_TEST_SUITE( graphics_toy )

#if 1
typedef GLuint idx_t;
#define GL_IDX_TYPE GL_UNSIGNED_INT
#else
typedef GLushort idx_t;
#define GL_IDX_TYPE GL_UNSIGNED_SHORT
#endif

UT_TEST_CASE( test1 )
{

  // initialize OpenGL
  avro_assert_msg( glfwInit() , "problem initializing OpenGL!" );
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
  //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // to make MacOS happy; should not be needed
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  // setup the primitives
  std::vector<GLfloat> coordinates = {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    2.0, 0.0, 0.0,

    0.0, 1.0, 0.0,
    1.0, 1.0, 0.0,
    2.0, 1.0, 0.0,

    0.0, 2.0, 0.0,
    1.0, 2.0, 0.0,
    2.0, 2.0, 0.0
  };

  idx_t indices[24] = {
    0,1,4,
    0,4,3,
    1,2,5,
    1,5,4,
    3,4,7,
    3,7,6,
    4,5,8,
    4,8,7
  };

  // create the window
  index_t width = 400;
  index_t height = width;
  GLFWwindow *window = glfwCreateWindow( width , height , "graphics toy" , NULL, NULL);
  if (!window) {
    const char* description;
    int code = glfwGetError(&description);

    if (description)
      printf("GLFW error (%d): %s\n",code,description);

    glfwTerminate();
    avro_assert_not_reached;
  }
  glfwMakeContextCurrent(window);

  // load GL functions
  gladLoadGL();
  dumpGLInfo();

  // define the basic shader
  ShaderProgram shader("basic");
  shader.use();

  // generate vertex array in which buffers will be saved
  GLuint vertex_array;
  GL_CALL( glGenVertexArrays( 1, &vertex_array ) );
  GL_CALL( glBindVertexArray(vertex_array) );

  // generate new buffers
  std::vector<GLuint> buffers(3);
  GL_CALL( glGenBuffers( buffers.size() , buffers.data() ) );
  GLuint position_buffer = buffers[0];
  GLuint triangle_buffer = buffers[1];

  // bind the position buffer
  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, position_buffer ) );
  GL_CALL( glBufferData(GL_ARRAY_BUFFER, coordinates.size() * sizeof(GLfloat), coordinates.data() , GL_STATIC_DRAW) );
  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, 0) );

  // bind the triangles
  GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, triangle_buffer ) );
  GL_CALL( glBufferData(GL_ELEMENT_ARRAY_BUFFER, 24 * sizeof(idx_t), indices , GL_STATIC_DRAW) );
  GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0) );

  // set up the view
  graphics::vec3 eye = {0,0,10};
  graphics::vec3 up = {0,1,0};
  graphics::vec3 center = {0.5,0.5,0.0};

  float fov = M_PI/4.0;
  float aspect = width/height;
  float znear  = 0.001;
  float zfar = 100;

  // define the transformation matrices, and set the MVP into the shader program
  graphics::mat4 perspective_matrix = glm::perspective( fov , aspect , znear , zfar );
  graphics::mat4 view_matrix = glm::lookAt( eye , center , up );
  graphics::mat4 model_matrix = glm::identity();
  graphics::mat4 mvp = perspective_matrix * view_matrix * model_matrix;
  shader.setUniform("MVP",mvp);

  // define which triangles we want to render
  std::vector<index_t> render_triangles = {0,2};

  // determine which indices need to be rendered
  std::vector<idx_t> render_indices(render_triangles.size()*3);
  index_t i = 0;
  for (index_t k = 0; k < render_triangles.size(); k++)
  for (index_t j = 0; j < 3; j++)
    render_indices[i++] = indices[ 3*render_triangles[k]+j ];
  GLuint start = *std::min_element( render_indices.begin() , render_indices.end() );
  GLuint end = *std::max_element( render_indices.begin() , render_indices.end() );
  GLuint count = render_indices.size();
  printf("start = %d, end = %d, count = %d\n",start,end,count);

  // bind the render triangles
  GLuint triangle_buffer_partial = buffers[2];
  GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, triangle_buffer_partial ) );
  GL_CALL( glBufferData(GL_ELEMENT_ARRAY_BUFFER, render_indices.size() * sizeof(idx_t), render_indices.data() , GL_STATIC_DRAW) );
  GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0) );

  // ensure we can capture the escape key being pressed
  glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
  glfwPollEvents();
  glfwSetCursorPos(window, width/2, height/2);

  glDisable(GL_CULL_FACE);
  glDisable(GL_DEPTH_TEST);

  GLuint active_buffer = triangle_buffer;

  // render loop
  while (true) {

    // grey-ish background color
    float col = 0.9;
    glClearColor (col,col,col, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    GL_CALL( glBindVertexArray(vertex_array) );
    GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, position_buffer  ) );
    GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) ); // enable attribute in position 0 which holds coordinates
    GL_CALL( glEnableVertexAttribArray(0) );
    GL_CALL( glBindBuffer( GL_ARRAY_BUFFER , 0 ) );

    // min index is 0, max index is 5 in first three triangles
    GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, active_buffer ) );
    //GL_CALL( glDrawRangeElements(GL_TRIANGLES, start , end , count , GL_UNSIGNED_INT , (void*)render_indices.data() ) );
    GL_CALL( glDrawElements(GL_TRIANGLES, 24 , GL_IDX_TYPE , 0 ) ) ;//(void*)indices ) );

    if (glfwWindowShouldClose(window)) break;
    if (glfwGetKey(window, GLFW_KEY_ESCAPE ) == GLFW_PRESS) break;
    if (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS) active_buffer = triangle_buffer;
    if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) active_buffer = triangle_buffer_partial;


    glfwSwapBuffers(window);
    glfwPollEvents();
  }

  glDeleteBuffers( buffers.size() , buffers.data() );
  glDeleteVertexArrays( 1 , &vertex_array );

  glfwTerminate();

}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(graphics_toy)
