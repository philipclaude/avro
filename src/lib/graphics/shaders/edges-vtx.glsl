#version 410
layout (location = 0 ) in vec3 a_Position;

out vec3 v_Position;

#if WITH_TESSELLATION == 0
out vec2 v_ParameterTess;
#endif

uniform int u_clip;
uniform vec3 u_clip_center;
uniform vec3 u_clip_normal;

flat out float v_clip;

void main()
{
  v_Position  = a_Position; // pass the position to the geometry shader
  gl_Position = vec4(a_Position,1.0);

  // determine if this vertex is visible
  v_clip = 1.0;
  if (u_clip == 2) {
    if (dot(a_Position - u_clip_center,u_clip_normal) < 0.0)
      v_clip = -1.0;
  }

  #if WITH_TESSELLATION == 0
  v_ParameterTess = vec2(0,0);
  #endif
}
