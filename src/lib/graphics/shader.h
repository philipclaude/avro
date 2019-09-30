#ifndef URSA_LIB_GRAPHICS_SHADER_H_
#define URSA_LIB_GRAPHICS_SHADER_H_

#include "graphics/gl.h"

#include "numerics/matrix.h"

#include <string>

namespace ursa
{

namespace numerics
{

typedef VectorS<2,float> vec2;
typedef VectorS<3,float> vec3;
typedef VectorS<4,float> vec4;
typedef MatrixS<3,3,float> mat3;
typedef MatrixS<4,4,float> mat4;

} // numerics

namespace graphics
{

enum GLSLShaderType
{
  VERTEX, FRAGMENT, GEOMETRY,
  TESS_CONTROL, TESS_EVALUATION
};

class ShaderProgram
{
public:
  ShaderProgram();

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

  void setUniform( const char *name, const numerics::vec2& v);
  void setUniform( const char *name, const numerics::vec3& v);
  void setUniform( const char *name, const numerics::vec4& v);
  void setUniform( const char *name, const numerics::mat3& m);
  void setUniform( const char *name, const numerics::mat4& m);

  void setUniform( const char *name, float val );
  void setUniform( const char *name, int val );
  void setUniform( const char *name, bool val );
  void setUniform( const char *name, int n, int *v);
  void setUniform( const char *name, int n, std::vector<int> v);

  void printActiveUniforms();
  void printActiveAttribs();

  int getUniformLocation(const char* name );
  bool fileExists( const std::string& fileName );

private:
  int handle_;
  bool linked_;
  std::string log_;

};

} // graphics

} // ursa

#endif
