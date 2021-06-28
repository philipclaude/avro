#include "unit_tester.hpp"

#include "graphics/gl.h"
#include "graphics/math.h"

using namespace avro;
using namespace avro::graphics;

UT_TEST_SUITE( graphics_toy )

UT_TEST_CASE( test1 )
{

  // initialize OpenGL
  avro_assert_msg( glfwInit() , "problem initializing OpenGL!" );
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
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

  GLuint indices[24] = {
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

  // define the basic shader
  ShaderProgram shader("basic");
  shader.use();

  // generate new buffers
  std::vector<GLuint> vbo(2);
  GL_CALL( glGenBuffers( vbo.size() , vbo.data() ) );
  GLuint vbo_position = vbo[0];
  GLuint vbo_triangles = vbo[1];

  // bind the triangles
  GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_triangles) );
  GL_CALL( glBufferData(GL_ELEMENT_ARRAY_BUFFER, 24 * sizeof(GLuint), indices , GL_STATIC_DRAW) );

  // bind the position buffer
  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, vbo_position) );
  GL_CALL( glBufferData(GL_ARRAY_BUFFER, coordinates.size() * sizeof(GLfloat), coordinates.data() , GL_STATIC_DRAW) );

  // bind buffers to the vertex attribute arrays
  GLuint vao_triangles;
  GL_CALL( glGenVertexArrays( 1, &vao_triangles ) );
  GL_CALL( glBindVertexArray(vao_triangles) );

  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, vbo_position ) );
  GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) ); // enable attribute in position 0 which holds coordinates
  GL_CALL( glEnableVertexAttribArray(0) );

  graphics::vec3 eye = {0,0,10};
  graphics::vec3 up = {0,1,0};
  graphics::vec3 center = {0.5,0.5,0.0};

  float fov = M_PI/4.0;
  float aspect = width/height;
  float znear  = 0.001;
  float zfar = 100;

  graphics::mat4 perspective_matrix = glm::perspective( fov , aspect , znear , zfar );
  graphics::mat4 view_matrix = glm::lookAt( eye , center , up );
  graphics::mat4 model_matrix = glm::identity();

  graphics::mat4 mvp = perspective_matrix * view_matrix * model_matrix;

  shader.setUniform("MVP",mvp);

  std::vector<index_t> render_triangles = {0,6,7};
  std::vector<GLuint> render_indices(render_triangles.size()*3);

  index_t i = 0;
  for (index_t k = 0; k < render_triangles.size(); k++)
  for (index_t j = 0; j < 3; j++)
    render_indices[i++] = indices[ 3*render_triangles[k]+j ];

  GLuint start = *std::min_element( render_indices.begin() , render_indices.end() );
  GLuint end = *std::max_element( render_indices.begin() , render_indices.end() );
  GLuint count = render_indices.size();
  printf("start = %u, end = %u, count = %lu\n",start,end,count);

  // render loop
  while (true) {

    // grey-ish background color
    float col = 0.9;
    glClearColor (col,col,col, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    GL_CALL( glBindVertexArray(vao_triangles) );
    
    // min index is 0, max index is 5 in first three triangles
    GL_CALL( glDrawRangeElements(GL_TRIANGLES, start , end , count , GL_UNSIGNED_INT , render_indices.data() ) );

    if (glfwWindowShouldClose(window)) break;
    glfwSwapBuffers(window);

    sleep(100);
    break;
  }

}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(graphics_toy)
