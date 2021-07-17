#version 410
layout (location = 0 ) in vec3 a_Position;

out vec3 v_Position;

#if WITH_TESSELLATION == 0
out vec2 v_ParameterTess;
#endif

void main()
{
  v_Position  = a_Position; // pass the position to the geometry shader
  gl_Position = vec4(a_Position,1.0);

  #if WITH_TESSELLATION == 0
  v_ParameterTess = vec2(0,0);
  #endif
}
