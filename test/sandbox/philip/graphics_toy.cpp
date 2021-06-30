#include "unit_tester.hpp"

#include "graphics/gl.h"
#include "graphics/math.h"

#include "graphics/colormap.h"

using namespace avro;
using namespace avro::graphics;

UT_TEST_SUITE( graphics_toy )

UT_TEST_CASE( test1 )
{
  // initialize OpenGL
  avro_assert_msg( glfwInit() , "problem initializing OpenGL!" );
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // to make MacOS happy; should not be needed
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

  std::vector<GLuint> indices = {
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
  std::vector<GLuint> buffers(5);
  GL_CALL( glGenBuffers( buffers.size() , buffers.data() ) );
  GLuint position_buffer = buffers[0];
  GLuint triangle_buffer = buffers[1];
  GLuint texture_buffer  = buffers[3];
  GLuint colormap_buffer = buffers[4];

  // bind the position buffer
  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, position_buffer ) );
  GL_CALL( glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * coordinates.size() , coordinates.data() , GL_STATIC_DRAW) );
  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, 0) );

  // bind the triangles
  GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, triangle_buffer ) );
  GL_CALL( glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * indices.size() , indices.data() , GL_STATIC_DRAW) );
  GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0) );

  // bind the solution values to the texture buffer
  std::vector<float> solution = {0,0,0,1,1,1};
  GL_CALL( glBindBuffer( GL_TEXTURE_BUFFER , texture_buffer) );
  GL_CALL( glBufferData( GL_TEXTURE_BUFFER , sizeof(float) * solution.size() , solution.data() , GL_STATIC_DRAW) );

  // generate a texture to hold the solution buffer
  GLuint texture;
  GL_CALL( glGenTextures( 1 , &texture) );
  GL_CALL( glActiveTexture( GL_TEXTURE0 + 0 ) );
  GL_CALL( glBindTexture( GL_TEXTURE_BUFFER , texture) );
  GL_CALL( glTexBuffer( GL_TEXTURE_BUFFER , GL_R32F , texture_buffer ) );

  // bind the colormap values to a buffer
  Colormap colormap;
  colormap.change_style("bwr");
  index_t ncolor = 256*3;
  GL_CALL( glBindBuffer( GL_TEXTURE_BUFFER , colormap_buffer) );
  GL_CALL( glBufferData( GL_TEXTURE_BUFFER , sizeof(float) * ncolor , colormap.data() , GL_STATIC_DRAW) );

  // generate a texture to hold the colormap buffer
  GLuint colormap_texture;
  GL_CALL( glGenTextures( 1 , &colormap_texture) );
  GL_CALL( glActiveTexture( GL_TEXTURE0 + 1 ) );
  GL_CALL( glBindTexture( GL_TEXTURE_BUFFER , colormap_texture) );
  GL_CALL( glTexBuffer( GL_TEXTURE_BUFFER , GL_R32F , colormap_buffer ) );

  // set up the view
  graphics::vec3 eye = {0,0,2};
  graphics::vec3 up = {0,1,0};
  graphics::vec3 center = {0.5,0.5,0.0};

  float fov    = M_PI/4.0;
  float aspect = width/height;
  float znear  = 0.001;
  float zfar   = 100;

  // define the transformation matrices, and set the MVP into the shader program
  graphics::mat4 perspective_matrix = glm::perspective( fov , aspect , znear , zfar );
  graphics::mat4 view_matrix = glm::lookAt( eye , center , up );
  graphics::mat4 model_matrix = glm::identity();
  graphics::mat4 mv = view_matrix * model_matrix;
  graphics::mat4 mvp = perspective_matrix * mv;
  shader.setUniform("u_ModelViewProjectionMatrix",mvp);

  // ensure we can capture the escape key being pressed
  glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
  glfwPollEvents();
  glfwSetCursorPos(window, width/2, height/2);

  glDisable(GL_CULL_FACE);
  glDisable(GL_DEPTH_TEST);

  // render loop
  while (true) {

    // grey-ish background color
    float col = 0.9;
    glClearColor (col,col,col, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // bind which attributes we want to draw
    GL_CALL( glBindVertexArray(vertex_array) );
    GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, position_buffer  ) );
    GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) ); // enable attribute in position 0 which holds coordinates
    GL_CALL( glEnableVertexAttribArray(0) );
    GL_CALL( glBindBuffer( GL_ARRAY_BUFFER , 0 ) );

    // bind the desired solution texture
    glActiveTexture(GL_TEXTURE0 + 0);
    GLint solution_location = glGetUniformLocation(shader.handle() , "solution");
    avro_assert( solution_location >= 0 );
    glUniform1i(solution_location, 0); // first sampler in fragment shader

    // bind the desired colormap
    glActiveTexture(GL_TEXTURE0 + 1);
    GLint colormap_location = glGetUniformLocation(shader.handle() , "colormap");
    avro_assert( colormap_location >= 0 );
    glUniform1i(colormap_location, 1); // second sampler in fragment shader

    // bind the buffer for the indices we want to draw
    GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, triangle_buffer ) );
    GL_CALL( glDrawElements(GL_TRIANGLES, 3 , GL_UNSIGNED_INT , 0 ) );

    // swap buffers and wait for user input
    glfwSwapBuffers(window);
    glfwPollEvents();

    // determine if we should exit the render loop
    if (glfwWindowShouldClose(window)) break;
    if (glfwGetKey(window, GLFW_KEY_ESCAPE ) == GLFW_PRESS) break;
  }

  // clean up
  glDeleteBuffers( buffers.size() , buffers.data() );
  glDeleteTextures( 1 , &colormap_texture );
  glDeleteTextures( 1 , &texture );
  glDeleteVertexArrays( 1 , &vertex_array );
  glfwTerminate();

}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(graphics_toy)
