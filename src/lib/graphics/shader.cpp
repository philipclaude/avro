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

const std::string __basic_vs_src__ =
#include "shaders/basic.vs"
;

const std::string __basic_fs_src__ =
#include "shaders/basic.fs"
;

const std::string __edge_fs_src__ =
#include "shaders/edge.fs"
;

const std::string __wv_vs_src__ =
#include "shaders/wv.vs"
;

const std::string __wv_fs_src__ =
#include "shaders/wv.fs"
;

namespace avro
{

namespace graphics
{

ShaderProgram::ShaderProgram( const std::string& name ) :
  handle_(-1),
  linked_(false),
  name_(name)
{
  if (name_=="basic")
  {
    avro_assert_msg( compile(name_.c_str(),__basic_vs_src__,__basic_fs_src__) , "error compiling basic shaders" );
  }
  else if (name_=="edge")
  {
    avro_assert_msg( compile(name_.c_str(),__wv_vs_src__,__wv_fs_src__) , "error compiling basic shaders" );
  }
  else if (name_=="wv")
  {
    avro_assert_msg( compile(name_.c_str(),__wv_vs_src__,__wv_fs_src__) , "error compiling basic shaders" );
  }
  else
    avro_implement;
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
ShaderProgram::setUniform( const char *name, mat4& m)
{
  if (!check()) return;
  int loc = getUniformLocation(name);
  if (loc >= 0)
  {
    GL_CALL( glUniformMatrix4fv(loc, 1, GL_FALSE, &m[0][0] ) );
  }
  else
    printf("uniform not set for program handle %d!!\n",handle_);
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

  const char* c_code = source.c_str();
  GL_CALL( glShaderSource( shaderHandle, 1, &c_code, NULL ) );

  // compile the shader
  GL_CALL( glCompileShader(shaderHandle ) );

  // check for errors
  int result;
  GL_CALL( glGetShaderiv( shaderHandle, GL_COMPILE_STATUS, &result ) );
  if (result == GL_FALSE)
  {
    // Compile failed, store log and return false
    int length = 0;
    log_ = "";
    GL_CALL( glGetShaderiv(shaderHandle, GL_INFO_LOG_LENGTH, &length ) );
    if (length > 0)
    {
      char* c_log = new char[length];
      int written = 0;
      GL_CALL( glGetShaderInfoLog(shaderHandle, length, &written, c_log) );
      log_ = c_log;
      delete [] c_log;
    }

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

void
ShaderLibrary::set_matrices( SceneGraph& scene )
{
  // go through all the active shaders and assign the MVP and normalMatrix
  std::map<std::string,ShaderProgram>::iterator it;
  for (it=shaders_.begin();it!=shaders_.end();it++)
  {
    ShaderProgram& shader = it->second;
    shader.use();
    shader.setUniform("MVP" , scene.mvp_matrix() );
    shader.setUniform("u_normalMatrix" , scene.normal_matrix() );
  }
}

void
ShaderLibrary::create()
{
  // create all the shaders
  shaders_.insert( { "wv" , ShaderProgram("wv") } );

  // set the uniforms
  ShaderProgram& shader = shaders_.at("wv");
  shader.use();
  shader.setUniform( "lightDir" , 0.0 , 0.3 , 1.0 );
  shader.setUniform( "wAmbient" , 0.25f );
  shader.setUniform( "xpar" , 1.0f );
  shader.setUniform( "conNormal" , 0. , 0. , 1. );
  shader.setUniform( "conColor" , 0. , 1. , 1. );
  shader.setUniform( "bacColor" , 0.5 , 0.5 , 0.5 );
  shader.setUniform( "wColor" , 1.0f );
  shader.setUniform( "bColor" , 0.0f );
  shader.setUniform( "wNormal" , 1.0f );
  shader.setUniform( "wLight" , 1.0f );
  shader.setUniform( "pointSize" , 2.0f );
  shader.setUniform( "picking" , 0 );
  shader.setUniform( "vbonum" , 0 );

  // TODO the rest of the shaders...
}

} // graphics

} // avro
