//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "common/error.h"

#include "graphics/shader.h"
#include "graphics/window.h"

#include "avro_config.h" // we need AVRO_SOURCE_DIR to find shaders

#include <fstream>

namespace avro
{

namespace graphics
{

#if AVRO_WITH_GL

std::string
get_shader_src( const std::string& filename ) {
  std::string content;
  std::ifstream file_stream(filename.c_str(), std::ios::in);

  if (!file_stream.is_open()) {
      std::cerr << "could not read file " << filename << ". file does not exist." << std::endl;
      return "";
  }

  std::string line = "";
  while (!file_stream.eof()) {
      std::getline(file_stream, line);
      content.append(line + "\n");
  }

  file_stream.close();
  return content;
}

ShaderProgram::ShaderProgram( const std::string& name , bool with_tess , const std::vector<std::string>& macros ) :
  handle_(-1),
  linked_(false),
  name_(name),
  macros_(macros),
  has_tessellation_shader_(with_tess)
{
  std::string base = AVRO_SOURCE_DIR + "/src/lib/graphics/shaders/"+ name;
  if (name_ == "wv") {
    std::string vtx_src = get_shader_src( AVRO_SOURCE_DIR + "/src/lib/graphics/shaders/wv-vtx.glsl" );
    std::string frg_src = get_shader_src( AVRO_SOURCE_DIR + "/src/lib/graphics/shaders/wv-frg.glsl" );
    avro_assert_msg( compile(name_.c_str(),vtx_src,frg_src) , "error compiling wv shader" );
  }
  else if (name == "points") {
    std::string vtx_src = get_shader_src( base + "-vtx.glsl" );
    std::string frg_src = get_shader_src( base + "-frg.glsl" );
    avro_assert_msg( compile(name_.c_str(),vtx_src,frg_src) , "error compiling points shader" );
  }
  else if (name == "basic" || name == "triangles" || name == "edges") {
    std::string vtx_src = get_shader_src( base + "-vtx.glsl" );
    std::string frg_src = get_shader_src( base + "-frg.glsl" );
    std::string geo_src = get_shader_src( base + "-geo.glsl" );

    if (with_tess) {
      std::string tcs_src = get_shader_src( base + "-tcs.glsl" );
      std::string tes_src = get_shader_src( base + "-tes.glsl" );
      avro_assert_msg( compile(name_.c_str(),vtx_src,frg_src,geo_src,tcs_src,tes_src) , "error compiling shader" );
    }
    else {
      avro_assert_msg( compile(name_.c_str(),vtx_src,frg_src,geo_src) , "error compiling basic shader" );
    }
  }
  else if (name == "particles") {
    std::string vtx_src = get_shader_src( base + "-vtx.glsl" );
    std::string frg_src = get_shader_src( base + "-frg.glsl" );
    avro_assert_msg( compile(name_.c_str(),vtx_src,frg_src) , "error compiling particle shader" );

    // i would like to use a geometry shader to render each particle as a sphere
    //std::string geo_src = get_shader_src( base + "-geo.glsl" );
    //avro_assert_msg( compile(name_.c_str(),vtx_src,frg_src,geo_src) , "error compiling particle shader" );

  }
  else {
    printf("unknown shader %s\n",name.c_str());
    avro_assert_not_reached;
  }
}

int
ShaderProgram::check()
{
  return true; // TODO
}

bool
ShaderProgram::link()
{
  if (!check()) return false;
  if (linked_) return true;
  if (handle_ <= 0) return false;

  GL_CALL( glLinkProgram(handle_) );

  int status = 0;
  GL_CALL( glGetProgramiv( handle_, GL_LINK_STATUS, &status) );
  if (status == GL_FALSE)
  {
    // Store log and return false
    int length = 0;
    log_ = "";

    GL_CALL( glGetProgramiv(handle_, GL_INFO_LOG_LENGTH, &length ) );

    if (length > 0)
    {
      char* c_log = new char[length];
      int written = 0;
      GL_CALL( glGetProgramInfoLog(handle_, length, &written, c_log) );
      log_ = c_log;
      delete [] c_log;
    }
    return false;
  }
  else
  {
    linked_ = true;
    return linked_;
  }
}

void
ShaderProgram::use()
{
  if (!check()) return;
  if (handle_ <= 0 || !linked_) return;
  //printf("handle = %d\n",handle_);
  GL_CALL( glUseProgram( handle_ ) );
}

std::string
ShaderProgram::log()
{
  return log_;
}

int
ShaderProgram::handle()
{
  return handle_;
}

bool
ShaderProgram::linked()
{
  return linked_;
}

void
ShaderProgram::bindAttribLocation( GLuint location, const char * name)
{
  if (!check()) return;
  GL_CALL( glBindAttribLocation(handle_, location, name) );
}

void
ShaderProgram::bindFragDataLocation( GLuint location, const char * name )
{
  if (!check()) return;
  GL_CALL( glBindFragDataLocation(handle_, location, name) );
}

void
ShaderProgram::setUniform( const char *name, float x, float y, float z)
{
  if (!check()) return;
  int loc = getUniformLocation(name);
  if (loc >= 0)
    GL_CALL( glUniform3f(loc,x,y,z) )
}

void
ShaderProgram::setUniform( const char *name, float x, float y, float z, float w)
{
  if (!check()) return;
  int loc = getUniformLocation(name);
  if (loc >= 0)
    GL_CALL( glUniform4f(loc,x,y,z,w) )
}

void
ShaderProgram::setUniform( const char *name, int n, float* v)
{
  float vf[10];
  if (!check()) return;
  int loc = getUniformLocation(name);
  if (loc >= 0)
  {
    if (n > 10) n = 10;
    for (int i=0; i<n;i++)
      vf[i] = v[i];
    GL_CALL( glUniform3fv(loc,n,vf) );
  }
}

void
ShaderProgram::setUniform( const char *name, const vec2& v)
{
  if (!check()) return;
  int loc = getUniformLocation(name);
  if (loc >= 0)
    GL_CALL( glUniform2f(loc,v[0],v[1]) )
}

void
ShaderProgram::setUniform( const char *name, const vec3& v)
{
  if (!check()) return;
  int loc = getUniformLocation(name);
  if (loc >= 0)
    GL_CALL( glUniform3f(loc,v[0],v[1],v[2]) )
}

void
ShaderProgram::setUniform( const char *name, const vec4& v)
{
  if (!check()) return;
  int loc = getUniformLocation(name);
  if (loc >= 0)
    GL_CALL( glUniform4f(loc,v[0],v[1],v[2],v[3]) )
}

void
ShaderProgram::setUniform( const char *name, mat3& m)
{
  if (!check()) return;
  int loc = getUniformLocation(name);
  if (loc >= 0)
    GL_CALL( glUniformMatrix3fv(loc, 1, GL_FALSE, &m[0][0] ) )
}

void
ShaderProgram::setUniform( const char *name, const mat4& m)
{
  if (!check()) return;
  int loc = getUniformLocation(name);
  if (loc >= 0)
  {
    GL_CALL( glUniformMatrix4fv(loc, 1, GL_FALSE, &m[0][0] ) );
  }
  else
    printf("uniform %s not set for program handle %d!!\n",name,handle_);
}

void
ShaderProgram::setUniform( const char *name, float val )
{
  if (!check()) return;
  int loc = getUniformLocation(name);
  if (loc >= 0)
    GL_CALL( glUniform1f(loc, val) )
}

void
ShaderProgram::setUniform( const char *name, int val )
{
  if (!check()) return;
  int loc = getUniformLocation(name);
  if (loc >= 0)
    GL_CALL( glUniform1i(loc, val) )
}

void
ShaderProgram::setUniform( const char *name, bool val )
{
  if (!check()) return;
  int loc = getUniformLocation(name);
  if (loc >= 0)
    GL_CALL( glUniform1i(loc, val) )
}

void
ShaderProgram::printActiveUniforms()
{
  GLint nUniforms, size, location, maxLen;
  GLchar* name;
  GLsizei written;
  GLenum type;

  if (!check()) return;

  GL_CALL( glGetProgramiv( handle_, GL_ACTIVE_UNIFORM_MAX_LENGTH, &maxLen) );
  GL_CALL( glGetProgramiv( handle_, GL_ACTIVE_UNIFORMS, &nUniforms) );

  name = (GLchar *) malloc( maxLen );

  printf("\nLocation | Name\n");
  printf("------------------------------------------------\n");
  for (int i=0;i<nUniforms;i++)
  {
    GL_CALL( glGetActiveUniform( handle_, i, maxLen, &written, &size, &type, name ) );
    location = glGetUniformLocation(handle_, name);
    printf(" %-8d | %s\n",location, name);
  }
  free(name);
}

void
ShaderProgram::printActiveAttribs()
{
  GLint written, size, location, maxLength, nAttribs;
  GLenum type;
  GLchar * name;

  if (!check()) return;

  glGetProgramiv(handle_, GL_ACTIVE_ATTRIBUTE_MAX_LENGTH, &maxLength);
  glGetProgramiv(handle_, GL_ACTIVE_ATTRIBUTES, &nAttribs);

  name = (GLchar *) malloc( maxLength );

  printf("\nIndex | Name\n");
  printf("------------------------------------------------\n");
  for (int i=0;i<nAttribs;i++)
  {
    GL_CALL( glGetActiveAttrib( handle_, i, maxLength, &written, &size, &type, name ) );
    location = glGetAttribLocation(handle_, name);
    printf(" %-5d | %s\n",location, name);
  }
  free(name);
}

bool
ShaderProgram::validate()
{
  if (!check()) return false;
  if (!linked()) return false;
  return true; // TODO

  GLint status;
  GL_CALL( glValidateProgram( handle_ ) );
  GL_CALL( glGetProgramiv( handle_, GL_VALIDATE_STATUS, &status ) );

  if (status == GL_FALSE)
  {
    // Store log and return false
    int length = 0;
    log_ = "";

    GL_CALL( glGetProgramiv(handle_, GL_INFO_LOG_LENGTH, &length ) );

    if (length > 0)
    {
      char* c_log = new char[length];
      int written = 0;
      GL_CALL( glGetProgramInfoLog(handle_, length, &written, c_log) );
      log_ = c_log;
      delete [] c_log;
    }

    return false;
  }
  else
  {
    return true;
  }
}

int
ShaderProgram::getUniformLocation(const char* name )
{
  if (!check()) return 0;
  return glGetUniformLocation(handle_, name);
}

bool
ShaderProgram::compileShaderFromString( const std::string& source, GLSLShaderType type )
{
  if (!check()) return false;

  if (handle_ <= 0)
  {
    handle_ = glCreateProgram();
    if (handle_ == 0)
    {
      log_ = "Unable to create shader program.\n";
      return false;
    }
  }

  GLuint shaderHandle = 0;

  switch (type)
  {
    case VERTEX:
      shaderHandle = glCreateShader(GL_VERTEX_SHADER);
      break;
    case FRAGMENT:
      shaderHandle = glCreateShader(GL_FRAGMENT_SHADER);
      break;
    case GEOMETRY:
      shaderHandle = glCreateShader(GL_GEOMETRY_SHADER);
      break;
    case TESS_CONTROL:
      shaderHandle = glCreateShader(GL_TESS_CONTROL_SHADER);
      break;
    case TESS_EVALUATION:
      shaderHandle = glCreateShader(GL_TESS_EVALUATION_SHADER);
      break;
    default:
      return false;
  }

  // add the macros after #version
  std::size_t idx = source.find("#version");
  idx = source.find("\n",idx);

  std::string source1 = source.substr(0,idx);
  std::string source2 = "\n";
  for (index_t k = 0; k < macros_.size(); k++)
    source2 += macros_[k] + "\n";
  std::string source3 = source.substr(idx+1,source.size());

  std::string final_source = source1 + source2 + source3;

  const char* c_code = final_source.c_str();
  GL_CALL( glShaderSource( shaderHandle, 1, &c_code, NULL ) );

  // compile the shader
  GL_CALL( glCompileShader(shaderHandle ) );

  // check for errors
  int result;
  GL_CALL( glGetShaderiv( shaderHandle, GL_COMPILE_STATUS, &result ) );
  if (result == GL_FALSE) {

    // compile failed, store log and return false
    int length = 0;
    log_ = "";
    GL_CALL( glGetShaderiv(shaderHandle, GL_INFO_LOG_LENGTH, &length ) );
    if (length > 0) {
      char* c_log = new char[length];
      int written = 0;
      GL_CALL( glGetShaderInfoLog(shaderHandle, length, &written, c_log) );
      log_ = c_log;
      delete [] c_log;
    }

    std::cout << final_source << std::endl;

    return false;
  }
  else
  {
    // compile succeeded, attach shader and return true
    GL_CALL( glAttachShader(handle_, shaderHandle) );
    return true;
  }
}

bool
ShaderProgram::compile(const char *name,
                     const std::string& vs,
                     const std::string& fs,
                     const std::string& gs,
                     const std::string& tcs,
                     const std::string& tes )
{
  if (!check()) return false;

  if (vs.empty() || fs.empty())
  {
   printf("GLSL programm error: %s: vertex and fragment shaders are mandatory to compile a GLSL programm \n",name);
   exit(1);
  }

  bool error = false;
  if (!compileShaderFromString(vs,VERTEX))
  {
     printf("GLSL programm error: %s: vertex shader(%s) failed to compile !\n%s",name,vs.c_str(),log().c_str());
     error = true;
  }
  if (!tcs.empty() && !tes.empty()  )
  {
    if (!compileShaderFromString(tcs,TESS_CONTROL))
    {
       printf("GLSL programm error: %s: tesselation control shader(%s) failed to compile !\n%s",name,tcs.c_str(),log().c_str());
       error = true;
    }
    if (!compileShaderFromString(tes,TESS_EVALUATION))
    {
       printf("GLSL programm error: %s: tesselation evaluation shader(%s) failed to compile !\n%s",name,tes.c_str(),log().c_str());
       error = true;
    }
  }
  if (!gs.empty())
  {
   if (!compileShaderFromString(gs,GEOMETRY))
   {
       printf("GLSL programm error: %s: geometry shader(%s) failed to compile !\n%s",name,gs.c_str(),log().c_str());
       error = true;
   }
  }
  if (!compileShaderFromString(fs,FRAGMENT))
  {
     printf("GLSL programm error: %s: fragment shader(%s) failed to compile !\n%s",name,fs.c_str(),log().c_str());
     error = true;
  }

  static const char* varyings[] = {"gl_Position","v_Color"};
  if (name_=="wv")
    GL_CALL( glTransformFeedbackVaryings(handle(),2,varyings,GL_INTERLEAVED_ATTRIBS) );

  if (!link())
  {
     printf("GLSL programm error: %s: shader failed to link!\n%s",name,log().c_str());
     error = true;
  }
  avro_assert(!error);
  return true;
}

#endif

} // graphics

} // avro
