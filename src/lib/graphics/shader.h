//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_GRAPHICS_SHADER_H_
#define avro_LIB_GRAPHICS_SHADER_H_

#include "graphics/gl.h"
#include "graphics/math.h"

#include <string>

namespace avro
{

namespace graphics
{

#if AVRO_WITH_GL

enum GLSLShaderType
{
  VERTEX, FRAGMENT, GEOMETRY,
  TESS_CONTROL, TESS_EVALUATION
};

class ShaderProgram
{
public:
  ShaderProgram( const std::string& name , bool with_tess=false , const std::vector<std::string>& macros = {} );

  bool compileShaderFromFile( const char* filename, GLSLShaderType type );
  bool compileShaderFromFile( const std::string& filename, GLSLShaderType type );
  bool compileShaderFromString( const std::string& source, GLSLShaderType type );

  bool compile(const char *name, const char *vs, const char *fs, const char *gs = NULL, const char *tcs = NULL, const char *tes = NULL);
  bool compile(const char *name, const std::string& vs, const std::string& fs, const std::string& gs = std::string(), const std::string& tcs  = std::string(), const std::string& tes =  std::string());

  bool link();
  bool validate();
  void use();
  int check();

  std::string log();
  int handle();
  bool linked();

  void bindAttribLocation( GLuint location, const char* name);
  void bindFragDataLocation( GLuint location, const char* name );

  void setUniform( const char *name, float x, float y, float z);
  void setUniform( const char *name, float x, float y, float z, float w);
  void setUniform( const char *name, int n, float* v);

  void setUniform( const char *name, const vec2& v);
  void setUniform( const char *name, const vec3& v);
  void setUniform( const char *name, const vec4& v);
  void setUniform( const char *name, mat3& m);
  void setUniform( const char *name, const mat4& m);

  void setUniform( const char *name, float val );
  void setUniform( const char *name, int val );
  void setUniform( const char *name, bool val );
  void setUniform( const char *name, int n, int *v);
  void setUniform( const char *name, int n, std::vector<int> v);

  void printActiveUniforms();
  void printActiveAttribs();

  int getUniformLocation(const char* name );
  bool fileExists( const std::string& filename );

  bool has_tessellation_shader() const { return has_tessellation_shader_; }

private:
  int handle_;
  bool linked_;
  std::string log_;
  std::string name_;
  std::vector<std::string> macros_;
  bool has_tessellation_shader_;
};

#endif

} // graphics

} // avro

#endif
