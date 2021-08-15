//#version 410
layout (location = 0 ) in vec3 a_Position;
layout (location = 1 ) in float a_Color;

uniform mat4 u_ModelViewProjectionMatrix;

out float v_Color;

void main()
{
  gl_Position  = u_ModelViewProjectionMatrix*vec4(a_Position,1.0);
  gl_PointSize = 10.0;

  v_Color = a_Color;
}
